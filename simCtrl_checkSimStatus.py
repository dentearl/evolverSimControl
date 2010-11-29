#!/usr/bin/env python
"""
simCtrl_checkTreeStatus.py
dent earl, dearl (a) soe ucsc edu
1 feb 2010

This script will look at a simCtrl output directory that is currently
being simulated, the newick tree that is directing the simulation, the
step size of the sim, and then finally the directory structure in order
to determine how much work has taken place (in steps) and how much remains
both in total steps and in approximate wall clock steps (assuming perfect
parallelism of the simulation tree).

And maybe one day it will produce a pretty picture about the whole thing.
"""
##############################
import copy
import glob
import math
import os
import sys
import time
import cPickle
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import simulation.lib.libSimTree as LST
from simulation.lib.libSimControl import stem

####################
# GRAPH SYMBOLS
symPreStat='O'
symStat='$'
symPreCycle='o'
symCycle='+'
symNone='-'
symTrans='#'

def usage():
    print 'USAGE: %s --help' %(sys.argv[0])
    print __doc__
    sys.exit(2)

class Branch:
    """I use 'Branch' objects to keep track of the
    longest branch and its name.
    """
    def __init__(self):
        self.path=0.0
        self.name=''

class Status:
    """I use a 'Status' object to keep track of all of the
    detailed results and current status in a single object. I find it easier
    to pass around. The items in the object are mostly defined
    in the unpackData() function.
    The variables dictionary is used to store whether or not certain
    variables have already been calculated (and thus do not need to be
    recalculated)
    """
    def __init__(self):
        self.name=''
        self.variables={}

def totalTreeStepsFinder(nt, ss):
    """totalTreeStepsFinder() takes a newickTree and a stepSize
    and finds the total number of steps in the entire tree.
    """
    numSteps=0
    if(nt == None):
        return int(numSteps)
    numSteps = math.ceil(nt.distance / ss)
    return int(numSteps + totalTreeStepsFinder(nt.right, ss) + totalTreeStepsFinder(nt.left, ss))

def totalTreeLengthFinder(nt):
    """totalTreeLengthFinder() takes a newickTree and
    finds the total length of the entire tree.
    """
    treeLength=0
    if(nt == None):
        return treeLength
    treeLength  = nt.distance
    return treeLength + totalTreeLengthFinder(nt.right) + totalTreeLengthFinder(nt.left)

def totalTreeDepthStepsFinder(nt, ss):
    """totalTreeDepthFinder() takes a newickTree and
    finds the depest leaf and returns the distance to the leaf
    """
    treeLength=0
    if(nt == None):
        return treeLength
    treeLength  = math.ceil(nt.distance / ss)
    return treeLength + max(totalTreeDepthStepsFinder(nt.right, ss), totalTreeDepthStepsFinder(nt.left, ss))

def completedTreeStepsFinder(nt, ss, cycles, comStepsDict, prgStepsDict):
    """completedTreeStepsFinder() takes a newickTree,
    a stepSize, and a list of the present cycle directories,
    and it returns the number of completed steps in the tree.
    """
    if(nt == None):
        return 0
    myCycles=cycles[:]
    runDirs=[]
    for i in range(0,len(myCycles)):
        myCycles[i]  = os.path.basename(myCycles[i])
        (head, tail) = os.path.split(cycles[i])
        runDirs.append(head)
    runDir = os.path.commonprefix(runDirs)
    cyclesDict = dict(zip(myCycles, [1]*len(myCycles)))
    while(nt.distance > 0):
        nt.distance = nt.distance - ss
        if nt.distance < 0:
            nt.distance=0
        name = LST.nameTree(nt)
        if name in cyclesDict and (name not in comStepsDict):
            if not os.path.exists(os.path.join(runDir, name,'cycleInfo.xml')):
                continue
            infotree=ET.parse(os.path.join(runDir, name,'cycleInfo.xml'))
            if infotree.find('timestamps'):
                ts=infotree.find('timestamps')
                if not ts.find('main'):
                    continue
                tsc=ts.find('main')
                pointDict={'cycleStep_1_cycleMain_1_start':1,
                           'cycleStep_1_cycleMain_1_end':2,
                           'cycleStep_2_cycleMain_2_start':3,
                           'cycleStep_2_cycleMain_2_end':4,
                           'cycleStep_3_cycleMain_3_start':5,
                           'cycleStep_3_cycleMain_3_end':6,
                           'cycleStep_4_cycleMain_4_start':7,
                           'cycleStep_4_cycleMain_4_end':8}
                for p in pointDict:
                    if tsc.find(p):
                        if name in prgStepsDict:
                            if prgStepsDict[name] < pointDict[p]:
                                prgStepsDict[name]=pointDict[p]
                        else:
                            prgStepsDict[name]=pointDict[p]
                if not ts.find('stats'):
                    continue
                tss=ts.find('stats')
                pointDict={'cycleStep_5_cycleStats_1_start':9,
                           'cycleStep_5_cycleStats_1_end':10,
                           'cycleStep_6_cycleStats_2_start':11,
                           'cycleStep_6_cycleStats_2_end':12,
                           'cycleStep_7_cycleStats_3_start':13,
                           'cycleStep_7_cycleStats_3_end':14,
                           'cycleStep_8_cycleStats_4_start':15,
                           'cycleStep_8_cycleStats_4_end':16}
                for p in pointDict:
                    if tss.find(p):
                        if prgStepsDict[name] < pointDict[p]:
                            prgStepsDict[name]=pointDict[p]
                toBeDeleted=[]
                for n in prgStepsDict:
                    if prgStepsDict[n] == 16:
                        comStepsDict[n]=1
                        toBeDeleted.append(n)
                for n in toBeDeleted:
                    del prgStepsDict[n]
    
    if(nt.distance <= 0):
        return max(len(comStepsDict),
                   completedTreeStepsFinder(nt.right, ss, cycles,
                                            comStepsDict, prgStepsDict) ,\
                   completedTreeStepsFinder(nt.left, ss, cycles, comStepsDict,
                                            prgStepsDict))
    else:
        return len(comStepsDict)

def longestRemainingBranchFinder(nt, ss, cycles):
    """longestRemainingBranchFinder() takes a newickTree,
    a stepSize, and a list of the present cycle directories,
    and returns a Branch() object with the longest path (in steps)
    and the name of the leaf-est branch.
    """
    b = Branch()
    path=0.0
    if(nt == None):
        return b
    while(nt.distance > 0):
        nt.distance = nt.distance - ss
        if nt.distance < 0:
            nt.distance=0
        path = path + 1
        b.path = b.path + 1
        name = LST.nameTree(nt)
        for i in cycles:
            if (i.split('/')[-1] == name):
                path=0.0
                del b
                b = Branch()
                b.name = name
    r = longestRemainingBranchFinder(nt.right, ss, cycles)
    l = longestRemainingBranchFinder(nt.left, ss, cycles)
#    print 'name: %s %f, right: %s %f left: %s %f\n' %(b.name, b.path, r.name, r.path, l.name, l.path)
    if  r.path >= l.path:
        b.path = b.path + r.path
        if r.name != '':
            b.name = r.name
    elif l.path > r.path:
        b.path = b.path + l.path
        if l.name != '':
            b.name = l.name
    return b

def elapsedTreeTimeFinder(nt, ss, simDir, cycles, elapsedTimesDict,
                          comStepsDict, parentNode='simulationInfo.xml'):
    """elapsedTreeTimeFinder() takes a newickTree, a stepSize, the location of the
    simulation directory, a list of all the present cycle directories, and the parentNode
    and returns the total amount of machine time spent on the tree including all branches.
    """
    etime=0.0
    myCycles=cycles[:]
    runDirs=[]
    for i in range(0,len(myCycles)):
        myCycles[i] = os.path.basename(myCycles[i])
        (head, tail) = os.path.split(cycles[i])
        runDirs.append(head)
    runDir = os.path.commonprefix(runDirs)
    cyclesDict = dict(zip(myCycles, [1]*len(myCycles)))
    if nt == None:
        return etime
    if (os.path.exists(os.path.join(simDir, parentNode))):
        if os.path.exists(os.path.join(simDir, parentNode)):
            infoTree=ET.parse(os.path.join(simDir, parentNode))
            root=infoTree.getroot()
            tObj=root.find('timestamps')
            if tObj.find('start') != None:
                if tObj.find('start').find('epochUTC') != None:
                    parentNodeTime = float(tObj.find('start').find('epochUTC').text)
    else:
        return etime
    running=True
    prgTimesDict = {}
    for c in cyclesDict:
        if c in comStepsDict:
            if c not in elapsedTimesDict:
                elapsedTimesDict[c] = getCycleRunTime(os.path.join(runDir,c,'cycleInfo.xml'))
        else:
            prgTimesDict[c] = getCycleRunTime(os.path.join(runDir,c,'cycleInfo.xml'))
        
    return (elapsedTimesDict, prgTimesDict)

def getCycleRunTime(file):
    time = 0.0
    if os.path.exists(file):
        infoTree=ET.parse(file)
        root=infoTree.getroot()
        tObj=root.find('timestamps')
        if tObj != None:
            if 'endEpochUTC' in tObj.attrib:
                time = ( float(tObj.attrib['endEpochUTC']) -
                         float(tObj.attrib['startEpochUTC']) )
    return time

def parentNodeStrFinder(target, newickStr, ss, previous='root'):
    """
    """
    nt = newickTreeParser(newickStr, 0.0)
    name = copy.copy(previous)
    prev = copy.copy(name)
    while(nt.distance > 0 ):
        nt.distance = nt.distance - ss
        if nt.distance < 0:
            nt.distance = 0
        prev = copy.copy(name)
        name = LST.nameTree(nt)
#        print 'name is now: '+name+', we seek: '+target+', name was: '+prev
        if target == name:
            return prev
    if nt.left != None:
        ltree = LST.nameTree(nt.left)
#        print 'going left: '+ltree+ ', previous: '+name+', tree: '+nt.left.iD+', target: '+target
        if nt.left.iD == target:
            return prev
        left  = parentNodeStrFinder(target, printBinaryTree(nt.left, 1), ss, name)
    else:
        left = ''
    if nt.right != None:
        rtree = LST.nameTree(nt.right)
        if nt.right.iD == target:
            return prev
#        print 'going right: '+rtree+', previous: '+name+', tree: '+nt.right.iD+', target: '+target
        right = parentNodeStrFinder(target, printBinaryTree(nt.right, 1), ss, name)
    else:
        right = ''
    if left != '':
        return left
    elif right != '':
        return right
    else:
        return ''
     
def howLongFinder(simDir, cycles, useNow=True):
    """howLongFinder() takes the simulation directory, a list of the
    present cycles, and a boolean and returns the length of time the
    simulation has either been running (useNow=True) or ran (useNow=False)
    """
    infoTree=ET.parse(os.path.join(simDir, 'simulationInfo.xml'))
    root=infoTree.getroot()
    ts=root.find('timestamps')
    st=ts.find('start')
    startTime=float(st.find('epochUTC').text)
    #startTime = os.path.getmtime(os.path.join(simDir, 'timeStamp_simulationStart'))
    if useNow:
        return time.time() - startTime
    endTime=0
    for name in cycles:
        if os.path.exists(os.path.join(simDir, name,'cycleInfo.xml')):
            infoTree=ET.parse(os.path.join(simDir, name,'cycleInfo.xml'))
            root=infoTree.getroot()
            tObj=root.find('timestamps')
            if tObj != None:
                if 'endEpochUTC' in tObj.attrib:
                    if endTime < float(tObj.attrib['endEpochUTC']):
                        endTime = float(tObj.attrib['endEpochUTC'])
    return endTime - startTime
        

def prettyTime(t):
    """ prettyTime(t): input t comes in as seconds. Output is a str
    with units attached (secs, mins, hrs, etc.)
    """
    if(t < 210):
        return '%.2f secs' %(t)
    t = t/60.0
    if(t < 121):
        return '%.2f mins' %(t)
    t = t/60.0
    if(t < 25):
        return '%.2f hrs' %(t)
    t = t/24.0
    if(t < 28):
        return '%.2f days' %(t)
    t = t/7
    return '%.2f weeks' %(t)

def directoriesOnly(list):
    """directoriesOnly() takes a list of items from a directory
    and purges anything from the list that is not a directory.
    """
    for i in list:
        if (not os.path.isdir(i)):
            list.remove(i)

def str2link(s, dir, title=''):
    if dir == '':
        return ''
    else:
        if title:
            return '<a href="'+dir+'/'+s+'/" title="'+title+'">'
        else:
            return '<a href="'+dir+'/'+s+'/">'

def prettyTitle(n, s):
    t = prettyTime(s)
    return 'Cycle %s took %s.' %(n, t)

def drawText(nt, ss, totalTreeDepth, rootName, scale=4, comStepsDict={}, prgStepsDict={}, isHtml=False, dir=''):
    """drawText() is in contrast to some other drawFORMAT() function that
    hasn't been written, but maybe one day will. drawText() draws the
    current state of the simulation as a tree using ASCII characters.
    """
    treeDepth=totalTreeDepthStepsFinder(nt, ss)
    depthFirstWalk(nt, stepSize=ss, scale=scale, rootName=rootName, comStepsDict=comStepsDict, prgStepsDict=prgStepsDict, isHtml=isHtml, dir=dir)
    drawScaleBar(treeDepth, scale, rootName, isHtml)

def depthFirstWalk(nt, stepSize=0.001, depth=0, branch='root', rootName='root', overlaps={}, scale=4, comStepsDict={}, prgStepsDict={}, isHtml=False, dir=''):
    """depthFirstWalk() depthFirstWalk. Walks a binaryTree object and writes
    out an ASCII representation of the tree.
    You need to know which branch you've descended from, in order
    to properly draw the vertical branchs between sister nodes.
    If you are the left branch and you are spawning a right branch,
    that right branch will need to add in a '|' at the correct position
    to make the tree appear connected, as in the 'All' level line for
    Steve or the 'S-Z-C' level line for Zack below:
            |###|###|###|###|###|###|###|###|+--|---| Jing
     ###|###| All
            |   |###|###|##+| Steve
            |###| S-Z-C
                |   |###|###|###|###|##+| Zack
                |###| Zack-Chris
                    |###|###|###|###|###|###|###|+--|---| Chris
    overlaps: depths where a pipe should be placed
    depth: current depth
    """
    if scale < 1:
        return
    if nt == None:
        return 0
    originalLength = nt.distance
    offset=' '*( len(rootName) - scale )
    for i in range(1, int(depth)+1):
        if str(i) in overlaps:
            offset += ' '*scale+' |'
        elif str(i-1) in overlaps:
            offset += ' '*scale
        else:
            offset += ' '*scale+' '
    if branch == 'root':
        depth=1
        if isHtml:
            offset = ' '+str2link(rootName, dir)+rootName+'</a>|'
        else:
            offset = ' %s|' % rootName
    else:
        offset+= '|'
    if isHtml:
        stringCap='</a>|'
    else:
        stringCap='|'
    for i in range(1,int(math.ceil(originalLength/stepSize))+1):
        # print symNone for steps not started, symStat for complete steps, symCycle for partially complete steps
        nt.distance = nt.distance - stepSize
        if nt.distance<0:
            nt.distance=0
        name = LST.nameTree(nt)
        if name in comStepsDict:
            offset=offset+str2link(name, dir, title=prettyTitle(name,comStepsDict[name]))+scale*symStat+stringCap
        elif name in prgStepsDict:
            if prgStepsDict[name]<9:
                # these are the + symbols, the main steps
                offset=offset+str2link(name, dir)+\
                        int(scale* int(prgStepsDict[name] - prgStepsDict[name]%2)/8) * symCycle +\
                        int(scale/4 * int(prgStepsDict[name])%2) * symPreCycle+\
                        (scale - int(scale*int(prgStepsDict[name] - prgStepsDict[name]%2)/8) -
                         int(scale/4 *int(prgStepsDict[name])%2) )*symNone+\
                         stringCap
            else:
                # these are the # symbols, the stats steps
                offset=offset+str2link(name, dir)+\
                        int(scale* int(prgStepsDict[name]- 8 - (prgStepsDict[name]%2))/8) * symStat +\
                        int(scale/4 * int(prgStepsDict[name])%2) * symPreStat+\
                        (scale - int(scale*int(prgStepsDict[name] - 8-  (prgStepsDict[name]%2))/8) -
                         int(scale/4 *int(prgStepsDict[name])%2) )*symCycle+\
                         stringCap
        else:
            offset=offset+scale*symNone+stringCap
    if not nt.right or not nt.left:
        print '%s %s' %(offset, nt.iD)
        return 1
    nextOverlapsR=dict.copy(overlaps)
    nextOverlapsL=dict.copy(overlaps)
    if branch=='left':
        # any time you double back, you're going to have to draw a connecting
        # branch ('|').
        nextOverlapsR[str(int(depth))]=1
    if branch=='right':
        nextOverlapsL[str(int(depth))]=1
    left = depthFirstWalk(nt.left, stepSize=stepSize, rootName=rootName, depth= depth+ (math.ceil(originalLength/stepSize)), branch='left',
               overlaps=nextOverlapsL, scale=scale, comStepsDict=comStepsDict, prgStepsDict=prgStepsDict, isHtml=isHtml, dir=dir)
    ##############################
    # FINALLY, print out the line and the name of the end cycle:
    if left:
        if nt.iD != rootName:
            print '%s %s' %(offset, nt.iD)
        else:
            print '%s' %( offset )
    right= depthFirstWalk(nt.right, stepSize=stepSize, rootName=rootName, depth=depth + (math.ceil(originalLength/stepSize)), branch='right',
               overlaps=nextOverlapsR, scale=scale, comStepsDict=comStepsDict, prgStepsDict=prgStepsDict, isHtml=isHtml, dir=dir) # +1 accounts for the '|' symbol
    if ( not left ) and right:
        if nt.iD != rootName:
            print '%s %s' %( offset, nt.iD )
        else:
            print '%s' %( offset )
    if left or right:
        return 1

def drawScaleBar(numSteps, scale, rootName, isHtml):
    scaleBar=' '* ( len( rootName ) - scale + 1)
    for i in range(int(numSteps+1), 0, -1):
        if i >= 100:
            if not i%10:
                scaleBar+=(int(scale/2)-1 )*' '+str(i)+(scale-int(scale/2)-3)*' '+'|'
            else:
                scaleBar+=scale*' '+'|'
        elif i >= 10:
            if not i%5:
                if scale > 2:
                    scaleBar+=(int(scale/2))*' '+str(i)+(scale-int(scale/2)-2)*' '+'|'
                elif scale > 1:
                    scaleBar+=str(i)+(scale-int(scale/2)-2)*' '+'|'
                else:
                    scaleBar+=' |'
            else:
                scaleBar+=scale*' '+'|'
        else:
            if not i%2:
                scaleBar+=(int(scale/2))*' '+str(i)+(scale-int(scale/2)-1)*' '+'|'
            else:
                scaleBar+=scale*' '+'|'
    print scaleBar
    print '%3s %20s %3s %20s %3s %20s\n%3s %20s %3s %20s %3s %20s' %(symNone, 'unfinished', symPreCycle, 'main step started',
                                                symCycle, 'main step complete', symPreStat,
                                                'stat step started', symStat, 'stat step complete',
                                                symTrans, 'transalign')
    if isHtml:
        print '</pre>'

def timeHandler(status, options):
    parentNode = parentNodeStrFinder(status.longBranchSteps.name, options.inputNewick, options.stepSize)
    if os.path.exists(os.path.join(options.dir, parentNode)):
        curCycleElapsedTime = time.time() - os.path.getmtime(os.path.join(options.dir, parentNode))
    else:
        curCycleElapsedTime = time.time() - os.path.getmtime(os.path.join(options.dir, status.longBranchSteps.name))
    (status.elapsedTreeTimeDict, prgTimeDict) = elapsedTreeTimeFinder(options.inputNewick, options.stepSize, options.dir,
                                                       status.cycleDirs, status.elapsedTreeTimeDict, status.completedTreeSteps_dict)
    elapsedTreeTime = sum(status.elapsedTreeTimeDict.values()) + sum(prgTimeDict.values())
    status.treeTime = prettyTime(elapsedTreeTime)
    newickTree = newickTreeParser(options.inputNewick, 0.0)
    
    #####
    # 
    if status.completedTreeSteps == status.totalTreeSteps:
        # the simulation is complete
        status.elapsedTime   = prettyTime(howLongFinder(options.dir, status.cycleDirs, useNow=False))
    else:
        # simulation is in progress
        status.longBranchSteps.path = status.longBranchSteps.path+0.5 # longBranchSteps fails to count the current running cycle
        status.elapsedTime          = prettyTime(howLongFinder(options.dir, status.cycleDirs))
    if status.completedTreeSteps and status.totalTreeSteps:
        # check to make sure these values are not 0
        status.aveBranchTime = prettyTime(elapsedTreeTime/status.completedTreeSteps)
        if status.completedTreeSteps != status.totalTreeSteps:
            if curCycleElapsedTime > (elapsedTreeTime/status.completedTreeSteps):
                # don't subtract off more than one cycle, it makes the estimates come out weird.
                rt = ((elapsedTreeTime/status.completedTreeSteps) * (status.longBranchSteps.path-0.5))
            else:
                rt = ((elapsedTreeTime/status.completedTreeSteps) * (status.longBranchSteps.path+0.5)) - curCycleElapsedTime
            #print ' ett/ctt = '+str(elapsedTreeTime/completedTreeSteps)+' lbs.p = '+str(longBranchSteps.path+0.5)+' curElTime: '+str(curCycleElapsedTime)
            if rt < 0:
                status.remainingTime = 'Soon'
                status.estTimeOfComp = 'Soon'
                status.estTotalRunLength = prettyTime(howLongFinder(options.dir, status.cycleDirs))
            else:
                status.estTotalRunLength = prettyTime(rt + howLongFinder(options.dir, status.cycleDirs))
                status.remainingTime = prettyTime(rt)
                status.estTimeOfComp = time.strftime("%a, %d %b %Y %H:%M:%S (UTC) ", time.gmtime(rt + time.time()))
        else:
            status.estTotalRunLength = prettyTime(howLongFinder(options.dir, status.cycleDirs, useNow=False))
            status.remainingTime = 'Done'
            status.estTimeOfComp = 'Done'
            status.variables['Done'] = True
    else:
        status.aveBranchTime = '--'
        status.remainingTime = '--'
        status.estTimeOfComp = '--'
        status.estTotalRunLength = '--'
    if status.longBranchSteps.path:
        status.workingCycleString = '(%s->%s)' %(parentNode, status.longBranchSteps.name)
    else:
        status.workingCycleString = ''

def initOptions(parser):
    parser.add_option('-d', '--dir',dest='dir',
                      help='Parent directory.')
    parser.add_option('-x', '--drawText',dest='drawText', action='store_true',default=False,
                      help='Prints an ASCII representation of the current tree status.')
    parser.add_option('-c', '--scale',dest='scale', action="store", default=4,
                      type ='int', help='Scale of the text phylogeny.')
    parser.add_option('-e', '--curCycles',dest='curCycles', action='store_true',default=False,
                      help='prints out the list of currently running cycles.')
    parser.add_option('-f', '--stepStats',dest='stepStats', action='store_true',default=False,
                      help='prints out the statistics for cycle steps.')
    parser.add_option('-g', '--cycleStem',dest='cycleStem', action='store_true',default=False,
                      help='prints out a stem and leaf plot for completed cycle runtimes, in seconds.')
    parser.add_option('--cycleStemHours',dest='cycleStemHours', action='store_true',default=False,
                      help='prints out a stem and leaf plot for completed cycle runtimes, in hours.')
    parser.add_option('-j', '--cycleList',dest='cycleList', action='store_true',default=False,
                      help='prints out a list of all completed cycle runtimes.')
    parser.add_option('-i', '--html',dest='isHtml', action='store_true',default=False,
                      help='prints output in an HTML friendly way, maybe for a cgi, hmmmmmmm?')
    parser.add_option('-t', '--htmlDir',dest='htmlDir', default='',
                      help='prefix for html links.')
    parser.add_option('--barimg',dest='barimg', default='',
                      help='URL to image file to be stretched to make barplots for html output.')

def checkOptions( options ):
    defaultStepSize = 0.001
    if not options.dir:
        options.dir=os.getcwd()
        options.isHtml=True
    if not os.path.isdir(options.dir):
        sys.stderr.write('%s: Error, directory "%s" is not a directory!\n' % (sys.argv[0], options.dir))
        usage()
    options.dir=os.path.abspath(options.dir)
    if not os.path.exists(os.path.join(options.dir, 'simulationInfo.xml')):
        sys.stderr.write('%s: Error, unable to locate simulationInfo.xml.\n' % sys.argv[0])
        usage()
    infoTree = ET.parse(os.path.join(options.dir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick = treeObj.text
    nt= newickTreeParser( options.inputNewick, 0.0 )
    if nt.distance==0:
        if nt.iD == None:
            nt.iD = 'root'
            options.inputNewick = printBinaryTree(nt,1)
    rootNameObj = infoTree.find('rootDir')
    options.rootName = os.path.basename( rootNameObj.text )
    treeObj = infoTree.find('saveParent')
    sp=treeObj.text
    options.saveParent = {'true':True, 'false':False}.get(sp.lower()) # and you thought python couldn't be cryptic?
    #                  # converts the strings 'true', 'false' to the Logical type python wants
    options.htmlDir = options.htmlDir.rstrip('/')
    
    # specify the stepSize for the simulation
    treeObj = infoTree.find('stepSize')
    stepSize=treeObj.text
    options.stepSize=float(stepSize)

    if options.isHtml:
        options.drawText=True
        pureHTMLStart()
        initHTML()

def pureHTMLStart():
    print 'Content-type: text/html'
    print ''

def initHTML():
    print '<!doctype html>'
    print '<html>'
    print '''
<head>
<title>Simulation Status</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
<style type="text/css">
pre, .code {
        padding: 10px 15px;
        margin: 5px 0 15px;
        border-left: 5px solid #CCCCCC;
        background: #FFFFFF; 
        font: 1em/1.5 "Courier News", monospace;
        color: #333333;
        overflow : auto;
}
hist { color: red;}
</style>
</head>
<body bgcolor="#FFFFFFFF">
'''

def finishHTML():
    print '''
</body>
</html>
'''

def csDictToTimeList(cs_dict, runDir):
    """csDictToTimeList() takes a dict of completed steps [cs_dict] and
    their common parent directory [runDir] and returns a list of all of
    their runtimes in seconds.
    """
    cycleTimes=[]
    for c in cs_dict:
        cycleTimes.append(cs_dict[c])
    return cycleTimes
def csDictToTimeDict(cs_dict, runDir):
    """csDictToTimeDict() takes a dict of completed steps [cs_dict] and
    their common parent directory [runDir] and returns a dict with cycle
    names as keys and their runtimes in seconds as values.
    """
    cycleTimes={}
    for c in cs_dict:
        if cs_dict[c] == 1:
            infotree=ET.parse(os.path.join(runDir, c,'cycleInfo.xml'))
            ts=infotree.find('timestamps')
            cycleTimes[c] = ( float(ts.attrib['endEpochUTC']) -
                              float(ts.attrib['startEpochUTC']) )
        else:
            cycleTimes[c] = cs_dict[c]
    return cycleTimes

def stepStatsBreakdown(runDir, cs_dict, isHtml, timeList, microTimeList,
                       transAlignTimeList, csAlreadyAdded, barimg):
    """stepStatsBreakdown() takes a dictionary of completed cycle steps and
    then finds out how long each step has taken and gives us some stats about
    the distribution of time spent in each step.
    cs_dict is the completed steps dict
    """
    timeDict={'cycleStep_1_cycleMain_1':timeList[0],
              'cycleStep_2_cycleMain_2':timeList[1],
              'cycleStep_3_cycleMain_3':timeList[2],
              'cycleStep_4_cycleMain_4':timeList[3],
              'cycleStep_5_cycleStats_1':timeList[4],
              'cycleStep_6_cycleStats_2':timeList[5],
              'cycleStep_7_cycleStats_3':timeList[6],
              'cycleStep_8_cycleStats_4':timeList[7]}
    mainSteps=['cycleStep_1_cycleMain_1', 'cycleStep_2_cycleMain_2',
               'cycleStep_3_cycleMain_3', 'cycleStep_4_cycleMain_4']
    statSteps=['cycleStep_5_cycleStats_1','cycleStep_6_cycleStats_2',
              'cycleStep_7_cycleStats_3', 'cycleStep_8_cycleStats_4']
    microSteps=['command_0', 'command_1', 'command_2',
                'command_3', 'command_4', 'command_5',
                'command_6']
    microStepsToNames={'command_0':'evolver_evo',
                       'command_1':'evolver_cvt',
                       'command_2':'trf_wrapper',
                       'command_3':'echo',
                       'command_4':'mv *.dat',
                       'command_5':'mv *.outseq.rev',
                       'command_6':'mv *.aln.rev'}
    microTimeDict={'evolver_evo':microTimeList[0],
                   'evolver_cvt':microTimeList[1],
                   'trf_wrapper':microTimeList[2],
                   'echo':microTimeList[3],
                   'mv *.dat':microTimeList[4],
                   'mv *.outseq.rev':microTimeList[5],
                   'mv *.aln.rev':microTimeList[6],}
    transalignSteps = ['command_0', 'command_1']
    transAlignToNames={'command_0':'transalign 1',
                       'command_1':'transalign 2'}
    transAlignTimeDict={'transalign 1':transAlignTimeList[0],
                        'transalign 2':transAlignTimeList[1]}
    cycleTimes=csDictToTimeList(cs_dict, runDir)
    for c in cs_dict:
        if c in csAlreadyAdded:
            # Since we're pickling our objects after every run, we should
            # only collect data on runtimes once and then after store it.
            continue        
        infotree=ET.parse(os.path.join(runDir, c,'cycleInfo.xml'))
        ts=infotree.find('timestamps')
        tsm=ts.find('main')
        for m in mainSteps:
            tObjStart = tsm.find(m+'_start')
            tObjEnd   = tsm.find(m+'_end')
            if tObjEnd and tObjStart:
                timeDict[m].append( float(tObjEnd.find('epochUTC').text) -
                                    float(tObjStart.find('epochUTC').text))
        tss=ts.find('stats')
        for s in statSteps:
            tObjStart = tss.find(s+'_start')
            tObjEnd   = tss.find(s+'_end')
            if tObjEnd and tObjStart:
                timeDict[s].append( float(tObjEnd.find('epochUTC').text) -
                                    float(tObjStart.find('epochUTC').text))
        chroms = glob.glob(os.path.join(runDir,c,'logs','micro.*.info.xml'))
        for chr in chroms:
            infotree = ET.parse(chr)
            ts  = infotree.find('timestamps')
            tsm = ts.find('micro')
            for m in microSteps:
                mObj = tsm.find(m)
                oObj = mObj.find('elapsed')
                elapsed = oObj.text
                microTimeDict[microStepsToNames[m]].append(float(elapsed))
        if os.path.exists(os.path.join(runDir,c,'logs','trans.info.xml')):
            infotree=ET.parse(os.path.join(runDir,c, 'logs','trans.info.xml'))
            ts  = infotree.find('timestamps')
            tsm = ts.find('micro')
            for t in transalignSteps:
                tsc = tsm.find(t)
                tse = tsc.find('elapsed')
                elapsed = tse.text
                transAlignTimeDict[transAlignToNames[t]].append(float(elapsed))
        csAlreadyAdded[c] = True
    step=1
    if isHtml:
        nStr = '<i>n</i>'
    else:
        nStr = 'n'
    maxMeanTime = 0
    for t in timeDict:
        if mean(timeDict[t]) > maxMeanTime:
            maxMeanTime = mean(timeDict[t])
    if maxMeanTime > 0:
        if isHtml:
            print '<h4>Step Stats</h4>'
        print 'time spent in each step (%s = %d)' %(nStr, len(cs_dict))
        if isHtml:
            print '<table cellpadding="5">'
            print '<tr><th>%s</th><th>%s</th><th>%s</th><th>%s</th><th>%s</th></tr>' %('step', 'ave (s)', 'pretty', 'sd', 'barplot')
        else:
            print '%13s %12s %8s %s' %('ave (s)', 'pretty', 'sd', 'barplot')
        for m in mainSteps:
            hist = '#'*int(25*(mean(timeDict[m])/maxMeanTime))
            if isHtml:
                hist = '&#9744;'*int(round(25*(mean(timeDict[m])/maxMeanTime))) # &#9744; &#183;
                sys.stdout.write( '<tr><td align="right" width=25>%d)</td>' % step )
                sys.stdout.write( '<td align="right" width=75>%8.2f</td>' % mean(timeDict[m]) )
                sys.stdout.write( '<td align="right" width=100>%s</td>' % prettyTime(mean(timeDict[m])) )
                if variance(timeDict[m]):
                    sdStr = '%8.2f' %math.sqrt(variance(timeDict[m]))
                else:
                    sdStr = '-'
                sys.stdout.write( '<td align="right" width=60>%s</td>' % sdStr )
                if barimg:
                    sys.stdout.write('<td align=left><img src="%s" height=10 width=%d></td></tr>\n' % (barimg, int(round(200*mean(timeDict[m])/maxMeanTime))))
                else:
                    sys.stdout.write( '<td><hist>%s</hist></td></tr>\n' % hist )
            else:
                hist = '#'*int(round(25*(mean(timeDict[m])/maxMeanTime)))
                print '%3d) %8.2f %12s %8.2f %s' %(step, mean(timeDict[m]), prettyTime(mean(timeDict[m])), math.sqrt(variance(timeDict[m])), hist)
            step += 1
        for s in statSteps:
            if isHtml:
                hist = '&#9744;'*int(round(25*(mean(timeDict[s])/maxMeanTime)))
                sys.stdout.write( '<tr><td align="right" width=25>%d)</td>' % step )
                sys.stdout.write( '<td align="right" width=75>%8.2f</td>' % mean(timeDict[s]) )
                sys.stdout.write( '<td align="right" width=100>%s</td>' % prettyTime(mean(timeDict[s])) )
                if variance(timeDict[s]):
                    sdStr = '%8.2f' %math.sqrt(variance(timeDict[s]))
                else:
                    sdStr = '-'
                sys.stdout.write( '<td align="right" width=60>%s</td>' % sdStr)
                if barimg:
                    sys.stdout.write('<td align=left><img src="%s" height=10 width=%d></td></tr>\n' % (barimg, int(round(200*mean(timeDict[s])/maxMeanTime))))
                else:
                    sys.stdout.write( '<td><hist>%s</hist></td></tr>\n' % hist )
            else:
                hist = '#'*int(round(25*(mean(timeDict[s])/maxMeanTime)))
                print '%3d) %8.2f %12s %8.2f %s' %(step, mean(timeDict[s]), prettyTime(mean(timeDict[s])), math.sqrt(variance(timeDict[s])), hist)
            step += 1
        
            
        if (isHtml) and ( maxMeanTime > 0.0 ):
            sys.stdout.write( '<tr><td>%s</td>' % 'total' )
            sys.stdout.write( '<td align="right">%8.2f</td>' % mean(cycleTimes) )
            sys.stdout.write( '<td align="right">%s</td>' % prettyTime(mean(cycleTimes)) )
            if variance(cycleTimes):
                sdStr = '%8.2f' %math.sqrt(variance(cycleTimes))
            else:
                sdStr = '-'
            sys.stdout.write( '<td align="right">%s</td>' % sdStr)
            sys.stdout.write( '<td></td></tr>' )
            print '</table>'
            
        print 'time spent in each chromosome (%s = %d)' %(nStr, len(microTimeList[0]))
        i = 0
        if isHtml:
            print '<table cellpadding="5">'
            print '<tr><th>%s</th><th>%s</th><th>%s</th><th>%s</th><th>%s</th></tr>' %('step', 'ave (s)', 'pretty', 'sd', 'barplot')
        else:
            print '%13s %12s %8s %s' %('ave (s)', 'pretty', 'sd', 'barplot')
        for m in microSteps:
            if isHtml:
                hist = '&#9744;'*int(round(25*(mean(microTimeDict[microStepsToNames[m]])/maxMeanTime)))
                sys.stdout.write( '<tr><td align="right" width=25>2.%d)</td>' % i )
                sys.stdout.write( '<td align="right" width=75>%8.2f</td>' % mean(microTimeDict[microStepsToNames[m]]) )
                sys.stdout.write( '<td align="right" width=100>%s</td>' % prettyTime(mean(microTimeDict[microStepsToNames[m]])) )
                if variance(microTimeDict[microStepsToNames[m]]):
                    sdStr = '%8.2f' % math.sqrt(variance(microTimeDict[microStepsToNames[m]]))
                else:
                    sdStr = '-'
                sys.stdout.write( '<td align="right" width=60>%s</td>' % sdStr)
                if barimg:
                    sys.stdout.write('<td align=left><img src="%s" height=10 width=%d></td></tr>\n' % (barimg, int(round(200*mean(microTimeDict[microStepsToNames[m]])/maxMeanTime))))
                else:
                    sys.stdout.write( '<td><hist>%s</hist></td></tr>\n' % hist )
            else:
                hist = '#'*int(round(25*(mean(microTimeDict[microStepsToNames[m]])/maxMeanTime)))
                print '2.%d) %8.2f %12s %8.2f %s' %(i, mean(microTimeDict[microStepsToNames[m]]),
                                                     prettyTime(mean(microTimeDict[microStepsToNames[m]])),
                                                     math.sqrt(variance(microTimeDict[microStepsToNames[m]])),
                                                     hist)
            i += 1
        if (isHtml) and ( maxMeanTime > 0.0 ):
            thisStep = 'cycleStep_2_cycleMain_2'
            sys.stdout.write( '<tr><td>%s</td>' % 'total' )
            sys.stdout.write( '<td align="right">%8.2f</td>' % mean(timeDict[thisStep]) )
            sys.stdout.write( '<td align="right">%s</td>' % prettyTime(mean(timeDict[thisStep])) )
            if variance(timeDict[thisStep]):
                sdStr = '%8.2f' % math.sqrt(variance(timeDict[thisStep]))
            else:
                sdStr = '-'
            sys.stdout.write( '<td align="right">%s</td>' % sdStr)
            sys.stdout.write( '<td></td></tr>' )
            print '</table>'
            
            
        print 'time spent in each transalign (%s ~ %d)' %(nStr, len(cs_dict))
        # it's approximate because the first cycles away from the root will not have transalign 1 times,
        # they'll have cp times which are very fast and as such are not recorded here. So the real number is
        # a little bit less than this n.
        i = 1
        if isHtml:
            print '<table cellpadding="5">'
            print '<tr><th>%s</th><th>%s</th><th>%s</th><th>%s</th><th>%s</th></tr>' %('step', 'ave (s)', 'pretty', 'sd', 'barplot')
        else:
            print '%13s %12s %8s %s' %('ave (s)', 'pretty', 'sd', 'barplot')
        for t in transAlignTimeDict:
            if isHtml:
                hist = '&#9744;'*int(round(25*(mean(transAlignTimeDict[t])/maxMeanTime)))
                sys.stdout.write( '<tr><td align="right" width=25>t %d)</td>' % i )
                sys.stdout.write( '<td align="right" width=75>%8.2f</td>'    % mean(transAlignTimeDict[t]) )
                sys.stdout.write( '<td align="right" width=100>%s</td>'       % prettyTime(mean(transAlignTimeDict[t])))
                if variance(transAlignTimeDict[t]):
                    sdStr = '%8.2f' % math.sqrt(variance(transAlignTimeDict[t]))
                else:
                    sdStr = '-'
                sys.stdout.write( '<td align="right" width=60>%s</td>'    % sdStr )
                if barimg:
                    sys.stdout.write('<td align=left><img src="%s" height=10 width=%d></td></tr>\n' % (barimg, int(round(200*mean(transAlignTimeDict[t])/maxMeanTime))))
                else:
                    sys.stdout.write( '<td><hist>%s</hist></td></tr>\n' % hist )
            else:
                hist = '#'*int(round(25*(mean(transAlignTimeDict[t])/maxMeanTime)))
                sys.stdout.write('t %d)'% i)
                sys.stdout.write(' %8.2f'% mean(transAlignTimeDict[t]))
                sys.stdout.write(' %12s' % prettyTime(mean(transAlignTimeDict[t])))
                if variance(transAlignTimeDict[t]):
                    sdStr = '%8.2f' % math.sqrt(variance(transAlignTimeDict[t]))
                else:
                    sdStr = '-'
                sys.stdout.write(' %s'% sdStr )
                sys.stdout.write(' %s\n' % hist )
            i += 1
        if (isHtml) and ( maxMeanTime > 0.0 ):
            print '</table>'
            

def findStalledCycles(runDir, cs_dict, prg_dict, isHtml, htmlDir=''):
    """findStalledCycles() calculates the average time cycles take and the
    standard deviation of that time and then it looks for cycles that are
    currently *in progress* and reports in progress cycles that are three
    SDs or more beyond the mean time.
    """
    if len(cs_dict) < 10:
        return
    cycleTimes=csDictToTimeList(cs_dict, runDir)
    meanTime = mean(cycleTimes)
    sdTime   = math.sqrt(variance(cycleTimes))
    toDelete_dict={}
    for p in prg_dict:
        if not os.path.exists(os.path.join(runDir, p ,'cycleInfo.xml')):
            toDelete_dict[p] = True
            continue
        infotree=ET.parse(os.path.join(runDir, p,'cycleInfo.xml'))
        ts=infotree.find('timestamps')
        tObjStart = float(ts.attrib['startEpochUTC'])
        if (time.time() - tObjStart) >= 3*sdTime + meanTime:
            if isHtml:
                print '<br><font color="red"><h2>WARNING</h2> cycle [ %s%s</a> ] has taken %s, which is > mean+3*SD (%s). Cycle may be stalled!</font><br>' %(str2link(p, htmlDir),p, prettyTime(time.time()-tObjStart), prettyTime(3*sdTime+meanTime))
            else:
                print '\nWARNING -  cycle [ %s%s ] has taken %s, which is > mean+3*SD (%s). Cycle may be stalled!\n' %(str2link(p, htmlDir),p, prettyTime(time.time()-tObjStart), prettyTime(3*sdTime+meanTime))
    for p in toDelete_dict:
        del(prg_dict[p]) # if this file doesn't exist, there's a good change the cycle was removed.

def listCurrentCycles(runDir, cs_dict, prg_dict, isHtml, htmlDir=''):
    """listCurretnCycles() lists all currently running cycles, their current
    step, current step runtime and current cycle runtime.
    """
    stepArray=[None,
               'cycleStep_1_cycleMain_1',
               'cycleStep_2_cycleMain_2',
               'cycleStep_3_cycleMain_3',
               'cycleStep_4_cycleMain_4',
               'cycleStep_5_cycleStats_1',
               'cycleStep_6_cycleStats_2',
               'cycleStep_7_cycleStats_3',
               'cycleStep_8_cycleStats_4']
    if len(prg_dict):
        if isHtml:
            print '<h4>Currently running cycles</h4>'
            print '<table cellpadding="5">'
            print '<tr><th align="right">%s</th><th>%s</th><th>%s</th><th>%s</th></tr>' %('Cycle', 'Total Time', 'Step', 'Step Time')
        else:
            print 'Currently running cycles:'
            print '%30s %15s %5s %15s' %('Cycle', 'Total Time', 'Step', 'Step Time')
    toDelete_dict={}
    keyList = prg_dict.keys()
    keyList.sort()
    for p in keyList:
        if not os.path.exists(os.path.join(runDir, p,'cycleInfo.xml')):
            toDelete_dict[p] = True
            continue
        infotree=ET.parse(os.path.join(runDir, p,'cycleInfo.xml'))
        ts=infotree.find('timestamps')
        tObjStart = float(ts.attrib['startEpochUTC'])
        runTime= (time.time() - tObjStart)
        stepNum=(prg_dict[p]+1)/2
        if stepNum<=4:
            tsm=ts.find('main')
        else:
            tsm=ts.find('stats')
        if prg_dict[p]%2:
            tObjStart=tsm.find(stepArray[stepNum]+'_start')
            tObjStart=float(tObjStart.find('epochUTC').text)
            stepTime = (time.time() - tObjStart)
        else:
            stepTime = 0.0
            stepNum  = stepNum + 1
        if isHtml:
            print '<tr><td>%s%s</a></td><td>%s</td><td align="center">%d</td><td>%s</td></tr>' %(str2link(p, htmlDir),p, prettyTime(runTime), stepNum, prettyTime(stepTime))
        else:
            print '%30s %15s %5d %15s' %(p, prettyTime(runTime), stepNum, prettyTime(stepTime))
    if isHtml and len(prg_dict):
        print '</table>'
    for p in toDelete_dict:
        del(prg_dict[p])

def printStem(aList, isHtml, inHours):
    """printStem() sends a list to the stem() function
    """
    timeUnits='seconds'
    if inHours:
        timeUnits='hours'
        for i in range(0, len(aList)):
            # times are in seconds.
            aList[i] = aList[i] / 60 / 60
    if isHtml:
        print '<h4>Distribution of cycle runtimes (%s)</h4>' %(timeUnits)
        print '<pre>'
    else:
        print 'Distribution of cycle runtimes (%s):' %(timeUnits)
    stem(aList)
    if isHtml:
        print '</pre>'
    
def printSortedCycleTimes(cs_dict, isHtml, htmlDir):
    """printSortedCycleTimes() prints out a list of completed cycles,
    sorted by the cycle runtimes.
    """
    items = cs_dict.items()
    returnItems = [[v[1], v[0]] for v in items]
    returnItems.sort(reverse=True)
    if isHtml:
        print '<h4>List of cycle runtimes</h4>'
        print '<table cellpadding="5">'
        print '<tr><th>Time (s)</th><th></th><th>Cycle</th></tr>'
    else:
        print 'List of cycle runtimes:'
    for k,v in returnItems:
        if isHtml:
            print '<tr><td align="right">%10.2f</td><td>%s</td><td>%s%s</a></td></tr>'%(float(k), prettyTime(k), str2link(v,htmlDir),v)
        else:
            print '%10.2f %10s %30s'%(float(k), prettyTime(k), v)
    if isHtml:
        print '</table>'
    
    
def mean(a):
    if len(a) > 0:
        return float(sum(a)/len(a))
    else:
        return 0.0

def variance(a):
    """variance() calculates a variance from a list, 'a'
    """
    meanA = mean(a)
    tot=0.0
    for i in a:
        tot += (float(i) - meanA)**2.0
    if len(a)-1 > 0:
        return (tot/(len(a)-1))
    else:
        return 0.0

def unpackData(file):
    """unpackData() opens up the pickle of the last run and pulls out
    all the relevant data.
    """
    if not os.path.exists(file):
        status=Status()
    else:
        FILE = open(file, 'r')
        status = cPickle.load(FILE)
        FILE.close()
        if 'isDone' in status.variables:
            status.variables['dontUpdate'] = True
    return status

def collectData(options, status):
    """collectNewData() takes in the options and a Status object
    and then goes about seeing what new data has been generated since
    the last run of checkTreeStatus.py
    """
    # The tree parsing and directory globing goes in here.
    status.cycleDirs                = glob.glob(os.path.join(options.dir,'*'))
    directoriesOnly(status.cycleDirs)
    if 'totalTreeLength' not in status.variables:
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        status.totalTreeLength          = totalTreeLengthFinder(newickTree)
        status.variables['totalTreeLength'] = True
    if 'totalTreeDepthSteps' not in status.variables:
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        status.totalTreeDepthSteps      = totalTreeDepthStepsFinder(newickTree, options.stepSize)
        status.variables['totalTreeDepthSteps'] = True
    if 'totalTreeSteps' not in status.variables:
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        status.totalTreeSteps           = totalTreeStepsFinder(newickTree, options.stepSize)
        status.variables['totalTreeSteps'] = True
    if 'completedTreeSteps_dict' not in status.variables:
        status.completedTreeSteps_dict  = {}
        status.inProgressTreeSteps_dict = {}
        status.variables['completedTreeSteps_dict'] = True
        status.variables['inProgressTreeSteps_dict'] = True
    if 'elapsedTreeTimeDict' not in status.variables:
        status.elapsedTreeTimeDict = {}
        status.variables['elapsedTreeTimeDict'] = True
    if 'timeList' not in status.variables:
        status.timeList=[[], [], [], [], [] ,[] ,[] ,[]]
        status.microTimeList=[[],[],[],[],[],[],[]]
        status.transAlignTimeList=[[],[]]
        status.csAlreadyAdded={}
        status.variables['timeList'] = True
    newickTree = newickTreeParser(options.inputNewick, 0.0)
    status.completedTreeSteps = completedTreeStepsFinder(newickTree, options.stepSize,
                                                         status.cycleDirs, status.completedTreeSteps_dict,
                                                         status.inProgressTreeSteps_dict)
    status.completedTreeSteps_dict = csDictToTimeDict(status.completedTreeSteps_dict, options.dir)
    #####
    # now we find the longest continuous branch yet to be simulated
    # longBranchSteps is a Branch() object
    newickTree = newickTreeParser(options.inputNewick, 0.0)
    status.longBranchSteps = longestRemainingBranchFinder(newickTree, options.stepSize, status.cycleDirs) 
    if not status.longBranchSteps.path:
        status.longBranchSteps.name = ''
        
    #####
    # Time Handling!
    timeHandler(status, options)

def packData(status, file):
    """packData() stores all of the Status in the appropriate pickle file.
    """
    FILE = open(file, 'w')
    cPickle.dump(status, FILE, 2) # 2 is the format protocol, 2=binary
    FILE.close()
    
def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    status=unpackData(os.path.join(options.dir, '.simulationStatus.pickle'))
    if 'isDone' not in status.variables:
        collectData(options, status)
    if options.isHtml:
        print '<table><tr><td>\n'
        elmDiv = '</td><td>'
        rowDiv = '<td></tr><tr><td>'
    else:
        elmDiv = ''
        rowDiv = ''
    info1 = 'Tot. tree len: %f (%d steps of %s)' % (status.totalTreeLength, status.totalTreeSteps,
                                                    str(options.stepSize).rstrip('0'))
    info2 = 'Longest remaining branch: %.1f %s' % (status.longBranchSteps.path, status.workingCycleString)
    print '%s%s%s%s' %( info1,elmDiv, info2, rowDiv )
    info1 = 'Tot. stps taken: %3d (%2.2f %% complete)' %(status.completedTreeSteps,
                                                         100*(float(status.completedTreeSteps) / float(status.totalTreeSteps)))
    info2 = 'Elapsed CPU time: %12s (ave: %12s/step)' %(status.treeTime, status.aveBranchTime)
    print '%s%s%s%s' %(info1, elmDiv, info2, rowDiv)
                                                                                              
    if options.isHtml:
        status.remainingTimeStr='<b>'+status.remainingTime+'</b>'
    else:
        status.remainingTimeStr=status.remainingTime
    info1 = 'ETtC: %12s' % status.remainingTimeStr
    info2 = 'EToC: %12s' % status.estTimeOfComp
    info3 = 'Elapsed wall-clock: %12s' % status.elapsedTime
    info4 = 'ETRL: %12s' % status.estTotalRunLength
    print '%s%s%s%s\n%s%s%s\n' % (info1, elmDiv, info2, rowDiv,
                                  info3, elmDiv, info4)
    if options.isHtml:
        print '</td></tr></table>'
        print '<br>'
    print 'Generated at %s' %(time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime()))
    if options.isHtml:
        print '<br>'
    #####
    # Draw the Tree!
    if(options.drawText):
        if options.isHtml:
            print '<pre>'
        nt = newickTreeParser(options.inputNewick, 0.0)
        drawText(nt, options.stepSize, status.totalTreeDepthSteps, options.rootName, scale=options.scale,
                 comStepsDict=status.completedTreeSteps_dict, prgStepsDict=status.inProgressTreeSteps_dict,
                 isHtml=options.isHtml, dir=options.htmlDir)
    
    findStalledCycles(options.dir, status.completedTreeSteps_dict, status.inProgressTreeSteps_dict, options.isHtml, options.htmlDir)
    if options.curCycles or options.isHtml:
        listCurrentCycles(options.dir, status.completedTreeSteps_dict, status.inProgressTreeSteps_dict, options.isHtml, options.htmlDir)
    if options.stepStats or options.isHtml:
        stepStatsBreakdown(options.dir, status.completedTreeSteps_dict, options.isHtml, status.timeList,
                           status.microTimeList, status.transAlignTimeList, status.csAlreadyAdded, options.barimg)
    if options.cycleStem or options.cycleStemHours or options.isHtml:
        printStem( csDictToTimeList(status.completedTreeSteps_dict, options.dir) , options.isHtml, options.cycleStemHours)
    if options.cycleList or options.isHtml:
        printSortedCycleTimes(status.completedTreeSteps_dict, options.isHtml, options.htmlDir)
    
    if options.isHtml:
        finishHTML()

    if 'dontUpdate' not in status.variables:
        packData(status, os.path.join(options.dir, '.simulationStatus.pickle'))

if __name__ == "__main__":
    main()
