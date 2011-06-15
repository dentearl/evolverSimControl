#!/usr / bin / env python
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
##################################################
# Copyright (C) 2009-2011 by
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##################################################
import copy
import cPickle
import evolverSimControl.lib.libSimControl as lsc
import glob
import math
from optparse import OptionParser
import os
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
import sys
import time
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat # exception handling for empty xml

####################
# GRAPH SYMBOLS
# these symbols show up in the legend and on the drawn tree
symbolDict = {'none' : '-',
              'cycleProg' : 'o',
              'cycleDone' : '+',
              'done' : '#',
              4 : '1',
              5 : '4',
              6 : '2',
              7 : '5',
              8 : '6',
              9 : '3',
              10 : '7'}
symPreStat  = 'O'
symStat     = '$'
symPreCycle = 'o'
symCycle    = '+'
symNone     = '-'
symTrans    = '#'

class Branch:
    """ We use Branch class to keep track of the
    longest branch and its name.
    """
    def __init__(self):
        self.name = ''
        self.longestPath = 0.0
        self.longestChild = ''

class Status:
    """ We use a 'Status' object to keep track of all of the
    detailed results and current simulation status. The items in 
    the object are mostly defined in the unpackData() function.
    The variables dictionary is used to store whether or not certain
    variables have already been calculated (and thus do not need to be
    recalculated)
    """
    def __init__(self):
        self.name = ''
        self.variables = {}

class Step:
    """ The Step class keeps track of a single step and all the timing
    information in that step.
    """
    def __init__(self):
        self.name        = '' # os.basename() of the directory
        self.complete    = False
        self.startTime   = -1
        self.endTime     = -1
        self.elapsedTime = -1
        self.timeDict    = {}

def initOptions(parser):
    parser.add_option('--dir', dest = 'simDir',
                      help='Parent directory.')
    parser.add_option('--drawText', dest = 'drawText', action = 'store_true',
                      default = False,
                      help = ('Prints an ASCII representation of the current '
                              'tree status. default=%default'))
    parser.add_option( '--curCycles', dest = 'curCycles', action = 'store_true',
                      default = False,
                      help = ('prints out the list of currently '
                              'running cycles. default=%default'))
    parser.add_option('--stats', dest = 'stepStats', action = 'store_true',
                      default = False,
                      help = ('prints out the statistics for cycle steps. default=%default'))
    parser.add_option( '--cycleStem', dest = 'cycleStem', action = 'store_true',
                      default = False,
                      help=('prints out a stem and leaf plot for completed '
                            'cycle runtimes, in seconds. default=%default'))
    parser.add_option('--cycleStemHours', dest = 'cycleStemHours', action = 'store_true',
                      default = False,
                      help=('prints out a stem and leaf plot for completed '
                            'cycle runtimes, in hours. default=%default'))
    parser.add_option('--cycleList', dest = 'cycleList', action = 'store_true',
                      default = False,
                      help = ('prints out a list of all completed cycle '
                              'runtimes. default=%default'))
    parser.add_option('--html', dest = 'isHtml', action = 'store_true',
                      default = False,
                      help=('prints output in an HTML friendly way, maybe '
                            'for a cgi, hmmmmmmm? default=%default'))
    parser.add_option('--htmlDir', dest = 'htmlDir', default = '',
                      help = 'prefix for html links. default=%default')
    parser.add_option('--barimg', dest = 'barimg', default = '',
                      help = ('URL to image file to be stretched '
                              'to make barplots for html output. default=%default'))

def checkOptions(options, parser):
    if options.simDir is None:
        options.simDir = os.getcwd()
        options.isHtml = True
    if not os.path.isdir(options.simDir):
        parser.error('directory "%s" is not a directory!' % options.simDir)
    options.simDir = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        parser.error('unable to locate simulationInfo.xml in %s.' % options.simDir)
    options.scale = 4
    infoTree = timeoutParse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick = treeObj.text
    nt = newickTreeParser(options.inputNewick, 0.0)
    if nt.distance == 0:
        if nt.iD is None:
            nt.iD = 'root'
            options.inputNewick = printBinaryTree(nt,1)
    rootNameObj = infoTree.find('rootName')
    options.rootName = rootNameObj.text
    
    options.htmlDir = options.htmlDir.rstrip('/')
    
    # specify the stepLength for the simulation
    treeObj = infoTree.find('stepLength')
    stepLength = treeObj.text
    options.stepLength = float(stepLength)

    if options.isHtml:
        options.drawText = True

def totalTreeStepsFinder(nt, sl):
    """totalTreeStepsFinder() takes a newickTree and a stepLength
    and finds the total number of steps in the entire tree.
    """
    numSteps = 0
    if nt is None:
        return int(numSteps)
    numSteps = math.ceil(nt.distance / sl)
    return int(numSteps + totalTreeStepsFinder(nt.right, sl) + totalTreeStepsFinder(nt.left, sl))

def totalTreeLengthFinder(nt):
    """totalTreeLengthFinder() takes a newickTree and
    finds the total length of the entire tree.
    """
    treeLength = 0
    if nt is None:
        return treeLength
    treeLength = nt.distance
    return treeLength + totalTreeLengthFinder(nt.right) + totalTreeLengthFinder(nt.left)

def totalTreeDepthStepsFinder(nt, sl):
    """totalTreeDepthFinder() takes a newickTree and
    finds the depest leaf and returns the distance to the leaf
    """
    treeLength = 0
    if nt is None:
        return treeLength
    treeLength = math.ceil(nt.distance / sl)
    return treeLength + max(totalTreeDepthStepsFinder(nt.right, sl), totalTreeDepthStepsFinder(nt.left, sl))

def simStepUpdater(nt, sl, cycleDirs, stepsDict, options):
    """ simStepUpdater goes through the tree and the stepsDict and updates the 
    steps in stepsDict that need updating.
    """
    stepList = tree2stepList(nt, sl)
    for s in stepList:
        if s not in stepsDict:
            stepsDict[s] = Step()
            stepsDict[s].name = s
    for s in stepsDict:
        if stepsDict[s].complete:
            continue
        updateTimingInfo(stepsDict[s], options)
    
    return numCompletedSteps(stepsDict), stepsDict

def numCompletedSteps(stepsDict):
    c = 0
    for s in stepsDict:
        if stepsDict[s].complete:
            c += 1
    return c

def updateTimingInfo(s, options):
    """ takes a Step() object and updates all the timing info associated
    with that sim step.
    """
    if not os.path.exists(os.path.join(options.simDir, s.name, 'xml')):
        return
    # summary.xml
    infotree = timeoutParse(os.path.join(options.simDir, s.name, 'xml', 'summary.xml'))
    if infotree is not None:
        ts = infotree.find('timestamps')
        if ts is not None:
            if 'startEpochUTC' in ts.attrib:
                s.startTime = float(ts.attrib['startEpochUTC'])
            if 'endEpochUTC' in ts.attrib:
                s.endTime = float(ts.attrib['endEpochUTC'])
                s.elapsedTime = s.endTime - s.startTime
                s.complete = True
    # cycle.xml stats.xml
    for t in ['Cycle', 'Stats']:
        infotree = timeoutParse(os.path.join(options.simDir, s.name, 'xml', '%s.xml' % t.lower()))
        if infotree is not None:
            ts = infotree.find('timestamps')
            if ts is not None:
                for elm in ['start', 'end']:
                    key = '%sEpochUTC' % elm
                    if key in ts.attrib:
                        s.timeDict[t + key] = float(ts.attrib[key])
                for i in xrange(1, 5):
                    for p in ['_start', '_end']:
                        key = '%sStep%d%s' % (t, i, p)
                        elm = ts.find(key)
                        if elm is None:
                            continue
                        elmUtc = elm.find('epochUTC')
                        if elmUtc is None:
                            continue
                        s.timeDict[key] = float(elmUtc.text)
    # trans.xml
    infotree = timeoutParse(os.path.join(options.simDir, s.name, 'xml', 'transalign.xml'))
    if infotree is not None:
        ts = infotree.find('timestamps')
        if ts is not None:
            for elm in ['start', 'end']:
                key = '%sEpochUTC' % elm
                if key in ts.attrib:
                    s.timeDict['Transalign' + key] = float(ts.attrib[key])
                key = 'TransalignStep1_%s' % elm
                elm = ts.find(key)
                if elm is None:
                    continue
                elmUtc = elm.find('epocUTC')
                if elmUtc is None:
                    continue
                s.timeDict[key] = float(elmUtc.text)

def timeoutParse(filename, timeout = 1.0, retry = 0.2):
    """ timeoutParse() takes an xml filename and attemps to parse the xml
    every [retry] seconds for a total of [timeout] seconds. We use this because
    the simulation must lock the xml to write data.
    """
    s = 0.0
    tree = None
    while s < timeout:
        try:
            tree = ET.parse(filename)
            break
        except (expat.ExpatError, IOError):
            time.sleep(retry)
            s += retry
    return tree

def longestRemainingBranchFinder(nt, sl, stepsDict):
    """longestRemainingBranchFinder() takes a newickTree,
    a stepLength, and a list of the present cycle directories,
    and returns a Branch() object with the longest path (in steps)
    and the name of the leaf-est branch.
    """
    b = Branch()
    path = 0.0
    if nt is None:
        return b
    while nt.distance > 0:
        nt.distance -= sl
        if nt.distance < 0:
            nt.distance = 0
        path += 1
        b.longestPath += 1
        name = lsc.nameTree(nt)
        if b.longestChild == '':
            b.longestChild = name
        if stepsDict[name].complete:
            path = 0.0
            del b
            b = Branch()
            b.name = name
    r = longestRemainingBranchFinder(nt.right, sl, stepsDict)
    l = longestRemainingBranchFinder(nt.left, sl, stepsDict)

    if  r.longestPath >= l.longestPath:
        b.longestPath += r.longestPath
        if r.name != '':
            b.name = r.name
            b.longestChild = r.longestChild
    else:
        b.longestPath += l.longestPath
        if l.name != '':
            b.name = l.name
            b.longestChild = l.longestChild
    # print('id: %s name: %s %f, right: %s %f left: %s %f\n' 
    #       % (nt.iD, b.name, b.longestPath, r.name, r.longestPath, l.name, l.longestPath))
    return b

def elapsedTreeTimeFinder(elapsedTimesDict, stepsDict):
    """elapsedTreeTimeFinder() takes a newickTree, a stepLength, the location of the
    simulation directory, a list of all the present cycle directories, and the parentNode
    and returns the total amount of machine time spent on the tree including all branches.
    """
    prgTimesDict = {}
    totalElapsedTime = 0.0
    for s in stepsDict:
        if stepsDict[s].startTime != -1:
            if stepsDict[s].endTime != -1:
                totalElapsedTime += stepsDict[s].elapsedTime
            else:
                totalElapsedTime += time.time() - stepsDict[s].startTime
    for c in stepsDict:
        if stepsDict[c].complete:
            if c not in elapsedTimesDict:
                elapsedTimesDict[c] = stepsDict[c].elapsedTime
        else:
            prgTimesDict[c] = time.time() - stepsDict[c].startTime
    return (elapsedTimesDict, prgTimesDict)

def howLongSimulationFinder(simDir, cycles, useNow = True):
    """howLongFinder() takes the simulation directory, a list of the
    present cycles, and a boolean and returns the length of time the
    simulation has either been running (useNow = True) or ran (useNow = False)
    """
    infoTree = timeoutParse(os.path.join(simDir, 'simulationInfo.xml'))
    root = infoTree.getroot()
    ts = root.find('timestamps')
    st = ts.find('start')
    startTime = float(st.find('epochUTC').text)
    if useNow:
        return time.time() - startTime
    if ts.find('end'):
        return float(ts.find('end').find('epochUTC').text) - startTime
    return time.time() - startTime
        

def prettyTime(t):
    """ prettyTime(t): input t comes in as seconds. Output is a str
    with units attached (secs, mins, hrs, etc.)
    """
    if t < 120:
        return '%.2f secs' % (t)
    t = t / 60.0
    if t < 120:
        return '%.2f mins' % (t)
    t = t / 60.0
    if t < 25:
        return '%.2f hrs' % (t)
    t = t / 24.0
    if t < 28:
        return '%.2f days' % (t)
    t = t / 7.0
    return '%.2f weeks' % (t)

def cycleDirectoriesOnly(aList, options):
    """cycleDirectoriesOnly() takes a list of items from a directory
    and purges anything from the list that is not a cycle directory.
    """
    out = []
    for i in aList:
        if (os.path.isdir(i) and os.path.exists(os.path.join(i, 'xml')) and 
            os.path.basename(i) != options.rootName):
            out.append(i)
    return out

def str2link(s, directory, title=''):
    if directory == '':
        return ''
    else:
        if title:
            return '<a href="%s/%s/" title="%s">' % (directory, s, title)
        else:
            return '<a href="%s/%s/">' % (directory, s)

def prettyTitle(n, s):
    t = prettyTime(s)
    return 'Cycle %s took %s.' % (n, t)

def drawText(nt, options, sl, totalTreeDepth, rootName, scale = 4, 
             stepsDict = {}, isHtml = False, directory = ''):
    """drawText() is in contrast to some other drawFORMAT() function that
    hasn't been written, but maybe one day will. drawText() draws the
    current state of the simulation as a tree using ASCII characters.
    """
    treeDepth = totalTreeDepthStepsFinder(nt, sl)
    depthFirstWalk(nt, options, stepLength = sl, scale = scale, rootName = rootName, 
                   stepsDict = stepsDict, isHtml = isHtml, directory = directory)
    drawScaleBar(treeDepth, scale, rootName, isHtml)
    drawLegend()

def depthFirstWalk(nt, options, stepLength = 0.001, depth = 0, branch='root', 
                   rootName='root', overlaps = {}, scale = 4, stepsDict = {}, 
                   isHtml = False, directory = ''):
    """depthFirstWalk() depthFirstWalk. Walks a binaryTree object and writes
    out an ASCII representation of the tree.
    You need to know which branch you've descended from, in order
    to properly draw the vertical branchs between sister nodes.
    If you are the left branch and you are spawning a right branch,
    that right branch will need to add in a '|' at the correct position
    to make the tree appear connected, as in the 'All' level line for
    Steve or the 'S-Z-C' level line for Zack below:
              |####|####|####|####|####|####|####|7777|++--|----| Jing
     ####|####| All
              |   |####|####|####| Steve
              |###| S-Z-C
                  |    |####|####|####|####|####| Zack
                  |####| Zack-Chris
                       |####|####|####|####|####|####|5555|+---|----| Chris
    + Legend -------------------------------------------------------------------------------------------------+
    |   -            unfinished   o         cycle running   +            cycle done                           |
    |   1   stat run trans -run   2   stat -run trans run   3  stat -run trans done   4  stat done trans -run |
    |   5      stat & trans run   6   stat run trans done   7   stat done trans run   #         step complete |
    +---------------------------------------------------------------------------------------------------------+
    overlaps: depths where a pipe should be placed
    depth: current depth offset
    """
    if nt is None:
        return 0
    originalBranchLength = nt.distance
    offset = ' ' * (len(rootName) - scale)
    for i in xrange(1, int(depth)+1):
        if str(i) in overlaps:
            offset += ' ' * scale + ' |'
        elif str(i - 1) in overlaps:
            offset += ' ' * scale
        else:
            offset += ' ' * scale + ' '
    if branch == 'root':
        depth = 1
        if isHtml:
            offset = ' ' + str2link(rootName, directory) + rootName + '</a>|'
        else:
            offset = ' %s|' % rootName
    else:
        offset += '|'
    if isHtml:
        stringCap = '</a>|'
    else:
        stringCap = '|'
    for i in xrange(0, int(math.ceil(originalBranchLength / stepLength))):
        # print symNone for steps not started, symStat for complete steps, 
        # symCycle for partially complete steps
        nt.distance -= stepLength
        if nt.distance < 0:
            nt.distance = 0
        name = lsc.nameTree(nt)
        if stepsDict[name].complete:
            # complete steps get filled in
            offset = (offset + str2link(name, directory, 
                                        title = prettyTitle(name, stepsDict[name].elapsedTime)) + 
                      scale * symbolDict['done'] + stringCap)
        else:
            if stepsDict[name].startTime == -1:
                offset = offset + scale * symbolDict['none'] + stringCap
            else:
                statsValue = getTernaryValue(stepsDict[name], 'stats', options)
                transValue = getTernaryValue(stepsDict[name], 'trans', options)
                symbValue = combineTernaryValues(statsValue, transValue)
                if symbValue == 3:
                    # neither the stat nor trans has started
                    offset = (offset + str2link(name, directory) 
                              + scale * symbolDict['cycleDone'] + stringCap)
                elif symbValue > 3:
                    offset = (offset + str2link(name, directory) 
                              + scale * symbolDict[symbValue] + stringCap)
            # if ('StatsStep4_end' in stepsDict[name].timeDict 
            #     and 'TransalignendEpochUTC' in stepsDict[name].timeDict):
            #     # these are the + symbols, the main steps
                
            #     offset = (offset + str2link(name, directory)
            #               + int(scale * int(prgStepsDict[name] - prgStepsDict[name] % 2) / 8)
            #               * symbolDict['cycleDone']
            #               + int(scale / 4 * int(prgStepsDict[name]) % 2) * symbolDict['cycleProg']
            #               + (scale - int(scale * int(prgStepsDict[name] - prgStepsDict[name] % 2) / 8)
            #                  - int(scale / 4 * int(prgStepsDict[name]) % 2)) * symbolDict['none']
            #               + stringCap)
            # else:
            #     # these are the stats steps
            #     offset = (offset + str2link(name, directory)
            #               + int(scale * int(prgStepsDict[name] - 8 - (prgStepsDict[name] % 2)) / 8) 
            #               * symStat + int(scale / 4 * int(prgStepsDict[name]) % 2) 
            #               * '*'
            #               + (scale - int(scale * int(prgStepsDict[name] - 8 - (prgStepsDict[name] % 2)) / 8) 
            #                  - int(scale / 4 * int(prgStepsDict[name]) % 2)) 
            #               * symbolDict['cycleDone'] + stringCap)
        # else:
        #     offset = offset + scale * symbolDict['none'] + stringCap
    if not nt.right or not nt.left:
        print '%s %s' % (offset, nt.iD)
        return True
    nextOverlapsR = dict.copy(overlaps)
    nextOverlapsL = dict.copy(overlaps)
    if branch == 'left':
        # any time you double back, you're going to have to draw a connecting
        # branch ('|').
        nextOverlapsR[str(int(depth))] = 1
    if branch == 'right':
        nextOverlapsL[str(int(depth))] = 1
    left = depthFirstWalk(nt.left, options, stepLength = stepLength, rootName = rootName, 
                          depth = depth + (math.ceil(originalBranchLength / stepLength)), 
                          branch = 'left', overlaps = nextOverlapsL, scale = scale, 
                          stepsDict = stepsDict, isHtml = isHtml, directory = directory)
    ##############################
    # FINALLY, print out the line and the name of the end cycle:
    if left:
        if nt.iD != rootName:
            print '%s %s' % (offset, nt.iD)
        else:
            print '%s' % offset
    right = depthFirstWalk(nt.right, options, stepLength = stepLength, rootName = rootName, 
                           depth = depth + (math.ceil(originalBranchLength / stepLength)), 
                           branch = 'right', overlaps = nextOverlapsR, scale = scale, 
                           stepsDict = stepsDict, isHtml = isHtml, directory = directory)
    if not left and right:
        if nt.iD != rootName:
            print '%s %s' % (offset, nt.iD)
        else:
            print '%s' % (offset)
    if left or right:
        return 1

def getTernaryValue(s, t, options):
    """ s is a Step() object and t is either 'stats' or 'trans'
    returns a 0,1,2 based on the status of the simulation step
    """
    if t not in ['stats', 'trans']:
        raise RuntimeError('Unanticipated value %s' % t)
    value = 0
    for i in xrange(0, 5):
        if '%s%sStep%d_start' % (t[0].upper(), t[1:], i) in s.timeDict:
            value = 1
            if i == 4:
                value = 2
    return value

def combineTernaryValues(statsValue, transValue):
    """ both values are in {0, 1, 2}
    and this function maps them to a unique space of
    {0, ..., 11}
    """
    return statsValue + (3 * transValue + 3)

def drawScaleBar(numSteps, scale, rootName, isHtml):
    scaleBar=' ' * (len(rootName) - scale + 1)
    for i in xrange(int(numSteps+1), 0, -1):
        if i >= 100:
            if not i % 10:
                scaleBar += (int(scale / 2) - 1) * ' ' + str(i) + (scale - int(scale / 2) - 3) * ' ' + '|'
            else:
                scaleBar += scale * ' ' + '|'
        elif i >= 10:
            if not i % 5:
                if scale > 2:
                    scaleBar += (int(scale / 2)) * ' ' + str(i) + (scale - int(scale / 2) - 2) * ' ' + '|'
                elif scale > 1:
                    scaleBar += str(i) + (scale - int(scale / 2) - 2) * ' ' + '|'
                else:
                    scaleBar +=' |'
            else:
                scaleBar += scale * ' ' + '|'
        else:
            if not i % 2:
                scaleBar += (int(scale / 2)) * ' ' + str(i) + (scale - int(scale / 2) - 1) * ' ' + '|'
            else:
                scaleBar += scale * ' ' + '|'
    print '+ Height %s+' % ('-' * len(scaleBar))
    print scaleBar

def drawLegend():
    print '+ Legend %s+' % ('-' * 97)
    print ('| %3s %21s %3s %21s %3s %21s%26s |\n'
           '| %3s %21s %3s %21s %3s %21s %3s %21s |\n'
           '| %3s %21s %3s %21s %3s %21s %3s %21s |'
           % (symbolDict['none'], 'unfinished', symbolDict['cycleProg'], 'cycle running',
              symbolDict['cycleDone'], 'cycle done', ' ',
              '1', 'stat run trans -run', '2', 'stat -run trans run',
              '3', 'stat -run trans done', '4', 'stat done trans -run',
              '5', 'stat & trans run',
              '6', 'stat done trans run', '7', 'stat run trans done',
              '#', 'step complete',
              ))
    print '+%s+' % ('-' * 105)

def timeHandler(status, options):
    if status.longBranchSteps.name == '':
        parentNode = ''
    else:
        parentNode = os.path.basename(lsc.getParentDir(os.path.join(options.simDir, 
                                                                    status.longBranchSteps.name)))

    if os.path.exists(os.path.join(options.simDir, status.longBranchSteps.name)):
        curCycleElapsedTime = time.time() - os.path.getmtime(os.path.join(options.simDir, 
                                                                          status.longBranchSteps.name))
    else:
        curCycleElapsedTime = time.time() - os.path.getmtime(os.path.join(options.simDir, 
                                                                          status.longBranchSteps.longestChild))
    (status.elapsedTreeTimeDict, prgTimeDict) = elapsedTreeTimeFinder(status.elapsedTreeTimeDict, 
                                                                      status.stepsDict)
    elapsedTreeTime = sum(status.elapsedTreeTimeDict.values()) + sum(prgTimeDict.values())
    status.treeTime = prettyTime(elapsedTreeTime)
    newickTree = newickTreeParser(options.inputNewick, 0.0)
    
    if status.numCompletedSteps == status.numTotalSteps:
        # the simulation is complete
        status.elapsedTime = prettyTime(howLongSimulationFinder(options.simDir, status.cycleDirs, useNow = False))
    else:
        # simulation is in progress
        # longBranchSteps counts the current running cycle so we subtract 0.5
        status.longBranchSteps.longestPath -= 0.5 
        status.elapsedTime = prettyTime(howLongSimulationFinder(options.simDir, status.cycleDirs))
    if status.numCompletedSteps and status.numTotalSteps:
        # check to make sure these values are not 0
        status.aveBranchTime = prettyTime(elapsedTreeTime / status.numCompletedSteps)
        if status.numCompletedSteps != status.numTotalSteps:
            if curCycleElapsedTime > (elapsedTreeTime / status.numCompletedSteps):
                # don't subtract off more than one cycle, it makes the estimates come out weird.
                # rt is remainingTime
                rt = ((elapsedTreeTime / status.numCompletedSteps) * 
                      status.longBranchSteps.longestPath)
            else:
                rt = ((elapsedTreeTime / status.numCompletedSteps) * 
                      status.longBranchSteps.longestPath) - curCycleElapsedTime
            if rt < 0:
                status.remainingTime = 'Soon'
                status.estTimeOfComp = 'Soon'
                status.estTotalRunLength = prettyTime(howLongSimulationFinder(options.simDir, status.cycleDirs))
            else:
                status.estTotalRunLength = prettyTime(rt + 
                                                      howLongSimulationFinder(options.simDir, status.cycleDirs))
                status.remainingTime = prettyTime(rt)
                status.estTimeOfComp = time.strftime('%a, %d %b %Y %H:%M:%S (UTC) ', 
                                                     time.gmtime(rt + time.time()))
        else:
            # simulation is complete
            status.estTotalRunLength = prettyTime(howLongSimulationFinder(options.simDir, 
                                                                          status.cycleDirs, useNow = False))
            status.remainingTime = 'Done'
            status.estTimeOfComp = 'Done'
            status.variables['Done'] = True
    else:
        # numComSteps or numPrgSteps are 0, indeterminate result
        status.aveBranchTime = '--'
        status.remainingTime = '--'
        status.estTimeOfComp = '--'
        status.estTotalRunLength = '--'
    if status.longBranchSteps.longestPath != '':
        status.workingCycleString = '(%s->%s)' % (status.longBranchSteps.name, 
                                                  status.longBranchSteps.longestChild)
    else:
        status.workingCycleString = ''

def pureHtmlStart():
    print 'Content-type: text / html'
    print ''

def initHtml():
    print '<!doctype html>'
    print '<html>'
    print '''
<head>
<title>Simulation Status</title>
<meta http-equiv="Content-Type" content="text / html; charset = iso-8859-1">
<style type="text / css">
pre, .code {
        padding: 10px 15px;
        margin: 5px 0 15px;
        border-left: 5px solid #CCCCCC;
        background: #FFFFFF; 
        font: 1em / 1.5 "Courier News", monospace;
        color: #333333;
        overflow : auto;
}
hist { color: red;}
</style>
</head>
<body bgcolor="#FFFFFFFF">
'''

def finishHtml():
    print '''
</body>
</html>
'''

def stepsDictToTimeList(stepsDict):
    """stepsDictToTimeList() takes a dict of steps [stepsDict] and
    their common parent directory [runDir] and returns a list of all of
    their runtimes in seconds.
    """
    cycleTimes = []
    for c in stepsDict:
        if stepsDict[c].complete:
            cycleTimes.append(stepsDict[c].elapsedTime)
    return cycleTimes
def stepsDictToTimeDict(stepsDict):
    """stepsDictToTimeDict() takes a dict of completed steps [cs_dict] and
    their common parent directory [runDir] and returns a dict with cycle
    names as keys and their runtimes in seconds as values.
    """
    cycleTimes = {}
    for c in stepsDict:
        if stepsDict[c].complete:
            cycleTimes[c] = stepsDict[c].elapsedTime
        else:
            cycleTimes[c] = time.time() - stepsDict[c].startTime
    return cycleTimes

def stepStatsBreakdown(options, status):
    """stepStatsBreakdown() takes a dictionary of completed cycle steps and
    then finds out how long each step has taken and gives us some stats about
    the distribution of time spent in each step.
    cs_dict is the completed steps dict
    """
    print 'stats step breakdown'

def findStalledCycles(runDir, stepsDict, isHtml, htmlDir = ''):
    """findStalledCycles() calculates the average time cycles take and the
    standard deviation of that time and then it looks for cycles that are
    currently *in progress* and reports in progress cycles that are three
    SDs or more beyond the mean time.
    """
    if len(stepsDict) < 10:
        # if there are fewer than 10 cycles complete then we cant get a decent
        # enough estimate of the mean time nor the stdev to issue warnings.
        return
    cycleTimes = stepsDictToTimeList(stepsDict)
    meanTime = mean(cycleTimes)
    if meanTime == 0:
        return
    sdTime = math.sqrt(variance(cycleTimes))
    for p in stepsDict:
        if stepsDict[p].complete:
            continue
        if stepsDict[p].startTime == -1:
            continue
        if (time.time() - stepsDict[p].startTime) >= 3 * sdTime + meanTime:
            if isHtml:
                print('<br><span style="color:red"><h2>WARNING</h2> cycle [%s%s</a>] '
                      'has taken %s, which is > mean + 3 * stdev (%s). Cycle may be '
                      'stalled!</span><br>' 
                      % (str2link(p, htmlDir), p,
                         prettyTime(time.time() - stepsDict[p].startTime), 
                         prettyTime(3 * sdTime + meanTime)))
            else:
                print('WARNING - cycle [%s%s] has taken %s, which is > mean '
                      '+ 3 * SD (%s). Cycle may be stalled!' 
                      % (str2link(p, htmlDir), p, 
                         prettyTime(time.time() - stepsDict[p].startTime), 
                         prettyTime(3 * sdTime + meanTime)))

def numRunningSteps(stepDict):
    """ takes a stepDict and returns the current 
    number of incomplete and running steps
    """
    n = 0
    for s in stepDict:
        if not stepDict[s].complete:
            if stepDict[s].startTime != -1:
                n += 1
    return n

def listCurrentCycles(runDir, stepDict, isHtml, options, htmlDir=''):
    """listCurretnCycles() lists all currently running cycles, their current
    step, current step runtime and current cycle runtime.
    """
    stepArray = [None]
    for c in ['Cycle', 'Stats']:
        for i in xrange(1,4):
            stepArray.append('%sStep%d_start' % (c, i))
            stepArray.append('%sStep%d_end' % (c, i))
    running = False
    if numRunningSteps(stepDict):
        running = True
        if isHtml:
            print '<h4>Currently running cycles</h4>'
            print '<table cellpadding="5">'
            print ('<tr><th align="right">%s</th><th>%s</th>'
                   '<th>%s</th><th>%s</th></tr>' % ('Cycle', 'Total Time', 'Step', 'Step Time'))
        else:
            print 'Currently running cycles:'
            print '%30s %15s %25s %15s' % ('Cycle', 'Total Time', 'Step', 'Step Time')
    for p in stepDict:
        if stepDict[p].complete:
            continue
        if stepDict[p].startTime == -1:
            continue
        runTime = time.time() - stepDict[p].startTime
        stepName, stepTime = currentStepInfo(stepDict[p], options.simDir)
        if isHtml:
            print('<tr><td>%s%s</a></td><td>%s</td><td align="center">%25s</td>'
                  '<td>%s</td></tr>' 
                  % (str2link(p, htmlDir), p, prettyTime(runTime), stepName, prettyTime(stepTime)))
        else:
            print '%30s %15s %25s %15s' % (p, prettyTime(runTime), stepName, prettyTime(stepTime))
    if isHtml and running:
        print '</table>'

def currentStepInfo(s, simDir):
    curStep = 'CycleStep'
    curTime = s.startTime
    for i in xrange(1, 5):
        key = 'CycleStep%d_start' % i
        if key in s.timeDict:
            curStep = 'CycleStep%d' % i
            curTime = time.time() - s.timeDict[key]
        else:
            break
    stats = False
    for i in xrange(1, 5):
        key = 'StatsStep%d_start' % i
        if key in s.timeDict:
            stats = True
            curStep = 'StatsStep%d' % i
            curTime = time.time() - s.timeDict[key]
        else:
            break
    if 'TransalignstartEpochUTC' in s.timeDict:
        transStep = 'transStep'
        if stats:
            curStep += ' ' + transStep
        else:
            curStep = transStep
        if 'TransalignendEpochUTC' not in s.timeDict:
            transTime = time.time() - s.timeDict['TransalignstartEpochUTC']
            if curTime < transTime:
                curTime = transTime
    return curStep, curTime

def printStem(aList, isHtml, inHours):
    """printStem() sends a list to the stem() function
    """
    timeUnits='seconds'
    if inHours:
        timeUnits='hours'
        for i in xrange(0, len(aList)):
            # times are in seconds.
            aList[i] = aList[i] / 60.0 / 60.0
    if isHtml:
        print '<h4>Distribution of cycle runtimes (%s)</h4>' % timeUnits
        print '<pre>'
    else:
        print 'Distribution of cycle runtimes (%s):' % timeUnits
    lsc.stem(aList)
    if isHtml:
        print '</pre>'
    
def printSortedStepTimes(stepsDict, isHtml, htmlDir):
    """printSortedStepTimes() prints out a list of completed cycles,
    sorted by the cycle runtimes.
    """
    completedTimesDict = {}
    for s in stepsDict:
        if stepsDict[s].complete:
            completedTimesDict[s] = stepsDict[s].elapsedTime
    items = completedTimesDict.items()
    returnItems = [[v[1], v[0]] for v in items]
    returnItems.sort(reverse = False)
    if isHtml:
        print '<h4>List of step runtimes</h4>'
        print '<table cellpadding="5">'
        print '<tr><th>Time (s)</th><th></th><th>Cycle</th></tr>'
    else:
        print 'List of step runtimes:'
        print '%12s %12s %30s' % ('Time (s)', '(Pretty)', 'Step name')
        print '%12s %12s %30s' % ('-' * 10, '-' * 10, '-' * 20)
    for k, v in returnItems:
        if isHtml:
            print('<tr><td align="right">%10.2f</td><td>%s</td><td>%s%s</a></td></tr>' 
                  % (float(k), prettyTime(k), str2link(v, htmlDir), v))
        else:
            print '%12.2f %12s %30s' % (float(k), prettyTime(k), v)
    if isHtml:
        print '</table>'
    
    
def mean(a):
    if len(a) > 0:
        return float(sum(a) / len(a))
    else:
        return 0.0

def variance(a):
    """variance() calculates a variance from a list, 'a'
    """
    meanA = mean(a)
    tot = 0.0
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
        status = Status()
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
    status.cycleDirs = glob.glob(os.path.join(options.simDir, '*'))
    status.cycleDirs = cycleDirectoriesOnly(status.cycleDirs, options)
    if 'totalTreeLength' not in status.variables:
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        status.totalTreeLength = totalTreeLengthFinder(newickTree)
        status.variables['totalTreeLength'] = True
    if 'totalTreeDepthSteps' not in status.variables:
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        status.totalTreeDepthSteps = totalTreeDepthStepsFinder(newickTree, options.stepLength)
        status.variables['totalTreeDepthSteps'] = True
    if 'numTotalSteps' not in status.variables:
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        status.numTotalSteps = totalTreeStepsFinder(newickTree, options.stepLength)
        status.variables['numTotalSteps'] = True
    if 'stepsDict' not in status.variables:
        status.stepsDict = {}
        status.variables['stepsDict'] = True
    if 'elapsedTreeTimeDict' not in status.variables:
        status.elapsedTreeTimeDict = {}
        status.variables['elapsedTreeTimeDict'] = True
    if 'timeDict' not in status.variables:
        status.timeDict = {}
        for c in ['Cycle', 'Stats']:
            for i in xrange(1,5):
                status.timeDict[c + 'Step' + str(i)] = []
        status.timeDict['TransalignStep1'] = []
        status.csAlreadyAdded = {}
        status.variables['timeDict'] = True
    newickTree = newickTreeParser(options.inputNewick, 0.0)
    newickTree.iD = options.rootName
    status.numCompletedSteps, status.stepsDict = simStepUpdater(newickTree, options.stepLength,
                                                                status.cycleDirs, status.stepsDict, 
                                                                options)
    #####
    # now we find the longest continuous branch yet to be simulated
    # longBranchSteps is a Branch() object
    newickTree = newickTreeParser(options.inputNewick, 0.0)
    newickTree.iD = options.rootName
    status.longBranchSteps = longestRemainingBranchFinder(newickTree, options.stepLength, status.stepsDict) 

    if not status.longBranchSteps.longestPath:
        status.longBranchSteps.name = ''
        
    #####
    # Time Handling!
    timeHandler(status, options)

def tree2stepList(nt, sl):
    """ takes a newick tree and a step length, returns
    a list of simulation step names
    """
    sList = []
    if nt is None:
        return sList
    while nt.distance > 0:
        nt.distance -= sl
        if nt.distance < 0:
            nt.distance = 0
        sList.append(lsc.nameTree(nt))
    return sList + tree2stepList(nt.right, sl) + tree2stepList(nt.left, sl)

def packData(status, filename):
    """packData() stores all of the Status in the appropriate pickle file.
    """
    f = open(filename, 'w')
    cPickle.dump(status, f) # 2 is the format protocol, 2 = binary
    f.close()
    
def printInfoTable(status, options):
    if options.isHtml:
        print '<table><tr><td>\n'
        elmDiv = '</td><td>'
        rowDiv = '<td></tr><tr><td>'
    else:
        elmDiv = ' '
        rowDiv = ' '
    info1 = ('Tot. tree len: %f (%d steps of %s)' 
             % (status.totalTreeLength, status.numTotalSteps,
                str(options.stepLength).rstrip('0')))
    info2 = ('Longest remaining branch: %.1f %s' 
             % (status.longBranchSteps.longestPath, status.workingCycleString))
    print '%s%s%s%s' % (info1,elmDiv, info2, rowDiv)
    info1 = ('Tot. stps taken: %d (%2.2f%% complete)' 
             % (status.numCompletedSteps,
                100 * (float(status.numCompletedSteps) / float(status.numTotalSteps))))
    info2 = '%sElapsed CPU time: %s (ave: %s / step)' % (' ' * 5, status.treeTime, status.aveBranchTime)
    print '%s%s%s%s' % (info1, elmDiv, info2, rowDiv)
                                                                                              
    if options.isHtml:
        status.remainingTimeStr = '<b>' + status.remainingTime + '</b>'
    else:
        status.remainingTimeStr = status.remainingTime
    info1 = 'ETtC: %31s ' % status.remainingTimeStr
    info2 = '%sEToC: %12s ' % (' ' * 4, status.estTimeOfComp)
    info3 = 'Elapsed wall-clock: %17s ' % status.elapsedTime
    info4 = '%sETRL: %12s' % (' ' * 4, status.estTotalRunLength)
    print('%s%s%s%s\n%s%s%s\n' 
          % (info1, elmDiv, info2, rowDiv, info3, elmDiv, info4))
    if options.isHtml:
        print '</td></tr></table>'
        print '<br>'
    print 'Generated at %s' % (time.strftime("%a, %d %b %Y %H:%M:%S (%Z)", time.localtime()))
    if options.isHtml:
        print '<br>'
    print 'Generated at %s' % (time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime()))
    if options.isHtml:
        print '<br>'
    else:
        print ' '

def printTree(status, options):
    #####
    # Draw the Tree!
    if options.drawText is not None:
        if options.isHtml:
            print '<pre>'
        nt = newickTreeParser(options.inputNewick, 0.0)
        drawText(nt, options, options.stepLength, status.totalTreeDepthSteps, options.rootName,
                  scale = options.scale, stepsDict = status.stepsDict,
                  isHtml = options.isHtml, directory = options.htmlDir)
        if options.isHtml:
            print '</pre>'
    
def printStats(status, options):
    findStalledCycles(options.simDir, status.stepsDict,
                       options.isHtml, options.htmlDir)
    if options.curCycles or options.isHtml:
        listCurrentCycles(options.simDir, status.stepsDict,
                           options.isHtml, options, options.htmlDir)
    if options.stepStats or options.isHtml:
        stepStatsBreakdown(options, status)
    if options.cycleStem or options.cycleStemHours or options.isHtml:
        printStem(stepsDictToTimeList(status.stepsDict),
                   options.isHtml, options.cycleStemHours)
    if options.cycleList or options.isHtml:
        printSortedStepTimes(status.stepsDict, options.isHtml, options.htmlDir)

def main():
    usage = ('usage: %prog --dir path/to/dir\n\n')
    parser = OptionParser(usage = usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, parser)

    status = unpackData(os.path.join(options.simDir, '.simulationStatus.pickle'))
    
    if 'isDone' not in status.variables:
        collectData(options, status)
        
    if options.isHtml:
        pureHtmlStart()
        initHtml()
        
    printInfoTable(status, options)
    printTree(status, options)
    printStats(status, options)
    
    if options.isHtml:
        finishHtml()

    if 'dontUpdate' not in status.variables:
        packData(status, os.path.join(options.simDir, '.simulationStatus.pickle'))

if __name__ == "__main__":
    main()
