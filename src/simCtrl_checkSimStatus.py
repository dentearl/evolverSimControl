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
import re
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
import subprocess
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

class SimNode:
    """ The SimNode object is used to store the binary tree of simulation
    nodes / directories. Useful for doing smart updating of timestamps.
    name should refer to a basename of a sim step directory inside of the 
    simulation dir. left right and parent are all dirs as well.
    """
    def __init__(self):
        self.name   = ''
        self.parent = None
        self.left   = None
        self.right  = None

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
    parser.add_option('--simDir', dest = 'simDir',
                      help = 'parent directory.')
    parser.add_option('--drawTree', '--drawText', dest = 'drawText', action = 'store_true',
                      default = False,
                      help = ('prints an ASCII representation of the current '
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
                      help = ('prints out a stem and leaf plot for completed '
                              'cycle runtimes, in seconds. default=%default'))
    parser.add_option('--cycleStemHours', dest = 'cycleStemHours', action = 'store_true',
                      default = False,
                      help = ('prints out a stem and leaf plot for completed '
                              'cycle runtimes, in hours. default=%default'))
    parser.add_option('--printChrTimes', dest = 'printChrTimes', action = 'store_true',
                      default = False,
                      help = ('prints a table of chromosome lengths (bp) and times '
                              '(sec) for intra chromosome evolution step (CycleStep2). '
                              'default=%default'))
    parser.add_option('--cycleList', dest = 'cycleList', action = 'store_true',
                      default = False,
                      help = ('prints out a list of all completed cycle '
                              'runtimes. default=%default'))
    parser.add_option('--html', dest = 'isHtml', action = 'store_true',
                      default = False,
                      help = ('prints output in HTML format for use as a cgi. '
                              'default=%default'))
    parser.add_option('--htmlDir', dest = 'htmlDir', default = '',
                      help = 'prefix for html links.')

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

def simStepUpdater(nt, sl, stepsDict, simNodeTree, status, options):
    """ simStepUpdater goes through the tree and the stepsDict and updates the 
    steps in stepsDict that need updating.
    """
    stepList = tree2stepList(nt, sl)
    for s in stepList:
        if s not in stepsDict:
            stepsDict[s] = Step()
            stepsDict[s].name = s
    simStepTreeUpdater(simNodeTree, stepsDict, status, options)
    return numCompletedSteps(stepsDict), stepsDict

def simStepTreeUpdater(snt, stepsDict, status, options):
    if snt is None:
        return
    started = False
    if snt.name != options.rootName:
        if snt.name not in stepsDict:
            raise RuntimeError('Unable to find expected SimNode snt.name = %s in stepsDict' % (snt.name))
        if not stepsDict[snt.name].complete:
            started = updateTimingInfo(stepsDict[snt.name], status, options)
        else:
            started = True
    else:
        started = True
    if started:
        # if the current node has not started, there is no 
        # sense in descending the tree
        simStepTreeUpdater(snt.left, stepsDict, status, options)
        simStepTreeUpdater(snt.right, stepsDict, status, options)

def numCompletedSteps(stepsDict):
    c = 0
    for s in stepsDict:
        if stepsDict[s].complete:
            c += 1
    return c

def updateTimingInfo(s, status, options):
    """ takes a Step() object and updates all the timing info associated
    with that sim step. If the Step has started running, returns True
    """
    updated = False
    if not os.path.exists(os.path.join(options.simDir, s.name, 'xml')):
        return updated
    # summary.xml
    infotree = timeoutParse(os.path.join(options.simDir, s.name, 'xml', 'summary.xml'))
    if infotree is not None:
        ts = infotree.find('timestamps')
        if ts is not None:
            if 'startEpochUTC' in ts.attrib:
                updated = True
                s.startTime = float(ts.attrib['startEpochUTC'])
            if 'endEpochUTC' in ts.attrib:
                updated = True
                s.endTime = float(ts.attrib['endEpochUTC'])
                s.elapsedTime = s.endTime - s.startTime
                s.complete = True
                updateChromLengthDict(s.name, status, options)
    # cycle.xml stats.xml
    for t in ['Cycle', 'Stats']:
        infotree = timeoutParse(os.path.join(options.simDir, s.name, 'xml', '%s.xml' % t.lower()))
        if infotree is not None:
            ts = infotree.find('timestamps')
            if ts is not None:
                for elm in ['start', 'end']:
                    key = '%sEpochUTC' % elm
                    if key in ts.attrib:
                        updated = True
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
                        updated = True
                        s.timeDict[key] = float(elmUtc.text)
    # chromosome run times
    chrs = glob.glob(os.path.join(options.simDir, s.name, 'xml', 'cycle.*.xml'))
    regex = r'cycle\.(.*)\.xml'
    pat = re.compile(regex)
    s.timeDict['chromosomes'] = {}
    for c in chrs:
        m = re.search(pat, os.path.basename(c))
        if m is None:
            raise RuntimeError('Regular expression %s failed to match chr xml file %s' % (regex, c))
        infotree = timeoutParse(c)
        if infotree is not None:
            ts = infotree.find('timestamps')
            if 'endEpochUTC' in ts.attrib:
                updated = True
                s.timeDict['chromosomes'][m.group(1)] = (float(ts.attrib['endEpochUTC']) - 
                                                         float(ts.attrib['startEpochUTC']))
    # trans.xml
    infotree = timeoutParse(os.path.join(options.simDir, s.name, 'xml', 'transalign.xml'))
    if infotree is not None:
        ts = infotree.find('timestamps')
        if ts is not None:
            for elm in ['start', 'end']:
                key = '%sEpochUTC' % elm
                if key in ts.attrib:
                    updated = True
                    s.timeDict['Transalign' + key] = float(ts.attrib[key])
                key = 'TransalignStep1_%s' % elm
                elm = ts.find(key)
                if elm is None:
                    continue
                elmUtc = elm.find('epochUTC')
                if elmUtc is None:
                    continue
                updated = True
                s.timeDict[key] = float(elmUtc.text)
    return updated

def updateChromLengthDict(name, status, options):
    """ calls evolver_cvt to get the length in bp for the chromosomes 
    in the 'name' step.
    should only be called once per simulation step, right when the step
    is recorded as 'complete'.
    """
    if 'chromosomeLengthsDict' not in status.variables:
        status.variables['chromosomeLengthsDict'] = True
        status.chromosomeLengthsDict = {}
    filename = os.path.join(options.simDir, name, 'seq.rev')
    cmd = [lsc.which('evolver_cvt')]
    cmd.append('-dumpchrids')
    cmd.append(filename)
    out = subprocess.Popen(cmd, stdout = subprocess.PIPE).communicate()[0]
    isChroms = False
    out = out.split('\n')
    for line in out:
        data = line.split()
        if data == []:
            continue
        if not isChroms:
            if len(data) != 4:
                continue
            else:
                if data[0].startswith('='):
                    isChroms = True
                    continue
        if isChroms:
            l = int(data[2])
            c = data[3]
            if c not in status.chromosomeLengthsDict:
                status.chromosomeLengthsDict[c] = {}
            if name not in status.chromosomeLengthsDict[c]:
                status.chromosomeLengthsDict[c][name] = l
            else:
                raise RuntimeError('sim step %s has identically named chromosomes "%s"' 
                                   % (name, c))

def timeoutParse(filename, timeout = 0.5, retry = 0.1):
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

def longestRemainingBranchFinder(snt, sl, stepsDict, options):
    """longestRemainingBranchFinder() takes a newickTree,
    a stepLength, and a list of the present cycle directories,
    and returns a Branch() object with the longest path (in steps)
    and the name of the leaf-est branch.
    """
    b = Branch()
    if snt is None:
        return b
    b.longestPath = 1.0
    b.longestChild = snt.name
    if snt.name != options.rootName:
        if stepsDict[snt.name].complete:
            path = 0.0
            del b
            b = Branch()
            b.name = snt.name
    else:
        b.name = options.rootName
        b.longestPath = 0.0
    l = longestRemainingBranchFinder(snt.left, sl, stepsDict, options)
    r = longestRemainingBranchFinder(snt.right, sl, stepsDict, options)
    if  l.longestPath >= r.longestPath:
        b.longestPath += l.longestPath
        if l.longestChild != '':
            b.longestChild = l.longestChild
        if l.name != '':
            b.name = l.name
    else:
        b.longestPath += r.longestPath
        if r.name != '':
            b.name = r.name
        if r.longestChild != '':
            b.longestChild = r.longestChild
    # print('id: %s name: %s %f, right: %s %f left: %s %f\n' 
    #       % (snt.name, b.name, b.longestPath, r.name, r.longestPath, l.name, l.longestPath))
    return b

def elapsedTreeTimeExtractor(elapsedTimesDict, stepsDict):
    """elapsedTreeTimeExtractor() takes a newickTree, a list of all the present cycle 
    directories, and the parentNode and returns the total amount of machine time 
    spent on the tree including all branches.
    """
    prgTimesDict = {}
    for c in stepsDict:
        if stepsDict[c].complete:
            if c not in elapsedTimesDict:
                elapsedTimesDict[c] = stepsDict[c].elapsedTime
        else:
            if stepsDict[c].startTime != -1:
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

def prettyLength(l):
    """ takes in an average base pairs (float) and kicks out a string
    """
    if not isinstance(l, float):
        raise RuntimeError('prettyLength is intended for floats, rewrite for %s' % l.__class__)
    if l > 10**6:
        return '%.2f Mb' % (l / float(10**6))
    if l > 10**3:
        return '%.2f Kb' % (l / float(10**6))
    return '%.2f' % l

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

def prettyTitle(n, s, done = True):
    t = prettyTime(s)
    if done:
        return 'Cycle %s took %s.' % (n, t)
    else:
        return 'Cycle %s has taken %s.' % (n, t)

def drawText(nt, options, sl, totalTreeDepth, scale = 4, 
             stepsDict = {}, isHtml = False, directory = ''):
    """drawText() is in contrast to some other drawFORMAT() function that
    hasn't been written, but maybe one day will. drawText() draws the
    current state of the simulation as a tree using ASCII characters.
    """
    treeDepth = totalTreeDepthStepsFinder(nt, sl)
    depthFirstWalk(nt, options, stepLength = sl, scale = scale, stepsDict = stepsDict, 
                   isHtml = isHtml, directory = directory)
    drawScaleBar(treeDepth, scale, options.rootName, isHtml)
    if isHtml:
        print '</pre>'
    drawLegend(isHtml)

def depthFirstWalk(nt, options, stepLength = 0.001, depth = 0, branch = 'root', 
                   overlaps = {}, scale = 4, stepsDict = {}, isHtml = False, 
                   directory = ''):
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
    offset = ' ' * (len(options.rootName) - scale)
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
            offset = ' ' + str2link(options.rootName, directory) + options.rootName + '</a>|'
        else:
            offset = ' %s|' % options.rootName
    else:
        offset += '|'
    if isHtml:
        linkCap = '</a>'
    else:
        linkCap = ''
    stringCap = '|'
    
    for i in xrange(0, int(math.ceil(originalBranchLength / stepLength))):
        nt.distance -= stepLength
        if nt.distance < 0:
            nt.distance = 0
        name = lsc.nameTree(nt)
        if stepsDict[name].complete:
            # complete steps get filled in
            link = str2link(name, directory, 
                            title = prettyTitle(name, stepsDict[name].elapsedTime))
            chunk = scale * symbolDict['done']
            assert len(chunk) == scale
            offset += link + chunk + linkCap + stringCap
        else:
            if stepsDict[name].startTime == -1:
                chunk = scale * symbolDict['none']
                assert len(chunk) == scale
                offset += chunk + stringCap
            else:
                statsValue = getTernaryValue(stepsDict[name], 'stats', options)
                transValue = getTernaryValue(stepsDict[name], 'trans', options)
                symbValue = combineTernaryValues(statsValue, transValue)
                if symbValue == 3:
                    # neither the stat nor trans has started
                    stepProg = getCycleStep(stepsDict[name])
                    lp = int(scale / 4 * stepProg)
                    if stepProg == 4:
                        stepDone = 4
                    else:
                        stepDone = stepProg - 1
                    ld = int(scale / 4 * stepDone)
                    link = str2link(name, directory, 
                                    title = prettyTitle(name, 
                                                        time.time() - stepsDict[name].startTime, 
                                                        done = False))
                    chunk = (ld * symbolDict['cycleDone']
                             + (lp - ld) * symbolDict['cycleProg']
                             + (scale - lp) * symbolDict['none'])
                    assert len(chunk) == scale
                    offset += link + chunk + linkCap + stringCap
                elif symbValue > 3:
                    link = str2link(name, directory, 
                                      title = prettyTitle(name, 
                                                          time.time() - stepsDict[name].startTime, 
                                                          done = False))
                    chunk = scale * symbolDict[symbValue]
                    assert len(chunk) == scale
                    offset += link + chunk + linkCap + stringCap
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
    left = depthFirstWalk(nt.left, options, stepLength = stepLength,
                          depth = depth + (math.ceil(originalBranchLength / stepLength)), 
                          branch = 'left', overlaps = nextOverlapsL, scale = scale, 
                          stepsDict = stepsDict, isHtml = isHtml, directory = directory)
    ##############################
    # FINALLY, print out the line and the name of the end cycle:
    if left:
        if nt.iD != options.rootName and nt.iD != None:
            print '%s %s' % (offset, nt.iD)
        else:
            print '%s' % offset
    right = depthFirstWalk(nt.right, options, stepLength = stepLength,
                           depth = depth + (math.ceil(originalBranchLength / stepLength)), 
                           branch = 'right', overlaps = nextOverlapsR, scale = scale, 
                           stepsDict = stepsDict, isHtml = isHtml, directory = directory)
    if not left and right:
        if nt.iD != options.rootName:
            print '%s %s' % (offset, nt.iD)
        else:
            print '%s' % (offset)
    if left or right:
        return 1

def getCycleStep(s):
    """ takes a Step() object and returns 1..4
    depending on which cycle step the Step is residing on
    """
    step = 1
    for i in xrange(1, 5):
        if 'CycleStep%d_start' % i in s.timeDict:
            step = i
    return step

def getTernaryValue(s, t, options):
    """ s is a Step() object and t is either 'stats' or 'trans'
    returns a 0,1,2 based on the status of the simulation step
    """
    if t not in ['stats', 'trans']:
        raise RuntimeError('Unanticipated value %s' % t)
    value = 0
    if t == 'stats':
        if 'StatsStep1_start' in s.timeDict:
            value = 1
    if t == 'trans':
        if 'TransalignstartEpochUTC' in s.timeDict:
            value = 1
    if t == 'stats' and 'StatsStep4_end' in s.timeDict:
        value = 2
    if t == 'transalign' and 'TransalignendEpochUTC' in s.timeDict:
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
    print '\n+ Height %s+' % ('-' * len(scaleBar))
    print scaleBar

def drawLegend(isHtml):
    if isHtml:
        print '<pre style="margin-left:2em;">'
    print '+ Legend %s+' % ('-' * 97)
    print ('| %3s %21s %3s %21s %3s %21s%26s |\n'
           '| %3s %21s %3s %21s %3s %21s %3s %21s |\n'
           '| %3s %21s %3s %21s %3s %21s %3s %21s |'
           % (symbolDict['none'], 'unstarted', symbolDict['cycleProg'], 'cycle running',
              symbolDict['cycleDone'], 'cycle done', ' ',
              '1', 'stat run trans -run', '2', 'stat -run trans run',
              '3', 'stat -run trans done', '4', 'stat done trans -run',
              '5', 'stat & trans run',
              '6', 'stat done trans run', '7', 'stat run trans done',
              '#', 'step complete',
              ))
    print '+%s+' % ('-' * 105)
    if isHtml:
        print '</pre>'

def timeHandler(status, options):
    curCycleElapsedTime = 0.0
    if status.longBranchSteps.name in status.stepsDict:
        if status.stepsDict[status.longBranchSteps.name].startTime != -1:
            curCycleElapsedTime = time.time() - status.stepsDict[status.longBranchSteps.name].startTime
    (status.elapsedTreeTimeDict, prgTimeDict) = elapsedTreeTimeExtractor(status.elapsedTreeTimeDict, 
                                                                         status.stepsDict)
    elapsedTreeTime = sum(status.elapsedTreeTimeDict.values()) + sum(prgTimeDict.values())
    status.elapsedTreeTimeStr = prettyTime(elapsedTreeTime)
    
    if status.numCompletedSteps == status.numTotalSteps:
        # the simulation is complete
        status.elapsedTime = howLongSimulationFinder(options.simDir, status.cycleDirs, 
                                                     useNow = False)
    else:
        # simulation is in progress
        status.elapsedTime = howLongSimulationFinder(options.simDir, status.cycleDirs)
    if status.numCompletedSteps and status.numTotalSteps:
        # check to make sure these values are not 0
        status.aveBranchTime = elapsedTreeTime / float(status.numCompletedSteps)
        status.aveBranchTimeStr = prettyTime(status.aveBranchTime)
        if status.numCompletedSteps != status.numTotalSteps:
            if curCycleElapsedTime > status.aveBranchTime:
                # don't subtract off more than one cycle, it makes the estimates come out weird.
                # rt is remainingTime
                rt = status.aveBranchTime * status.longBranchSteps.longestPath
            else:
                rt = status.aveBranchTime * status.longBranchSteps.longestPath - curCycleElapsedTime
            if rt < 0:
                status.remainingTimeStr = 'Soon'
                status.estTimeOfCompStr = 'Soon'
                status.estTotalRunLength = prettyTime(status.elapsedTime)
            else:
                status.estTotalRunLength = prettyTime(rt + status.elapsedTime)
                status.remainingTimeStr = prettyTime(rt)
                status.estTimeOfCompStr = time.strftime('%a, %d %b %Y %H:%M:%S (UTC) ', 
                                                     time.gmtime(rt + time.time()))
        else:
            # simulation is complete
            status.variables['isDone'] = True
            status.estTotalRunLength = prettyTime(status.elapsedTime)
            status.remainingTimeStr = 'Done'
            status.estTimeOfCompStr = 'Done'
    else:
        # numComSteps or numPrgSteps are 0, indeterminate result
        status.aveBranchTimeStr = '--'
        status.remainingTimeStr = '--'
        status.estTimeOfCompStr = '--'
        status.estTotalRunLength = '--'
    if status.longBranchSteps.longestPath != '':
        status.workingCycleString = '(%s --> %s)' % (status.longBranchSteps.name, 
                                                     status.longBranchSteps.longestChild)
    else:
        status.workingCycleString = ''

def pureHtmlStart():
    print 'Content-type: text / html'
    print ''

def initHtml(status):
    print '<!doctype html>'
    print '<html>'
    if 'isDone' in status.variables:
        refresh = ''
    else:
        # refresh = '<meta http-equiv="refresh" content="30;">'
        refresh = '''
<script type="text/javascript">
  function StartTime(){
    setTimeout("RefreshPage()", 30000);
  }
  function RefreshPage(){
    if(document.Reload.checkboxReload.checked){
      document.location.href = document.location.href;
    }
  }
  window.onload = StartTime();
</script>
'''
    print '''
<head>
<title>Simulation Status</title>
<meta http-equiv="Content-Type" content="text / html; charset = iso-8859-1">
%s

<style type="text/css">
pre, .code {
  padding: 10px 15px;
  margin: 5px 0 15px;
  border-left: 5px solid #CCCCCC;
  background: #FFFFFF; 
  font: 1em / 1.5 "Courier News", monospace;
  color: #333333;
  overflow : auto;
}
span.chrname {
  color:blue;
}
span.chrname:hover{
  background-color:#EFE1A1;
}
</style>
</head>
<body bgcolor="#FFFFFFFF">
''' % refresh

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
    if options.isHtml:
        print '<h3>Step Stats</h3>'
        print '<div style="margin-left:2em;">'
        print '<table border="1" bordercolor="#cccccc"><thead>'
        print ('<tr><th>%s</th><th>%s</th><th>%s</th></tr>' 
               % ('Cycles', 'Stats', 'Transalign'))
        print '</thead><tbody><tr>'
        for s in ['cycle', 'stats', 'trans']:
            print '<td valign="top">'
            printStatsSection(s, options, status)
            print '</td>'
        print '</tr></tbody></table></div>'
    else:
        print 'Step Stats breakdown'
        for s in ['cycle', 'stats', 'trans']:
            printStatsSection(s, options, status)

def printStatsSection(s, options, status):
    if s not in ['cycle', 'stats', 'trans']:
        raise TypeError('Unanticipated input to printStats: %s' % s)
    if s == 'cycle':
        printCycleStats(options, status)
    elif s == 'stats':
        printCycleStats(options, status, 'Stats')
    elif s == 'trans':
        printCycleStats(options, status, 'Transalign', 1)

def printCycleStats(options, status, pre = 'Cycle', length = 4):
    subList = []
    subTimesDict = {}
    for i in xrange(1, length + 1):
        subList.append('%sStep%d' % (pre, i))
        subTimesDict['%sStep%d' % (pre, i)] = []
    times = []
    for s in status.stepsDict:
        if '%sendEpochUTC' % pre in status.stepsDict[s].timeDict:
            times.append(status.stepsDict[s].timeDict['%sendEpochUTC' % pre] -
                         status.stepsDict[s].timeDict['%sstartEpochUTC' % pre])
        for u in subList:
            if '%s_end' % u in status.stepsDict[s].timeDict:
                subTimesDict[u].append(status.stepsDict[s].timeDict['%s_end' % u] - 
                                       status.stepsDict[s].timeDict['%s_start' % u])
    
    if options.isHtml:
        if pre == 'Cycle':
            extraCol = '<td style="text-align:center">--</td>'
            extraColHeader = '<th>&mu; Len.</th>'
        else:
            extraCol = ''
            extraColHeader = ''
        print('<table cellpadding="5"><thead><tr><th>Name</th><th>'
              '<span style="font-style:italic">n</span></th>%s'
              '<th>&mu; time (s)</th><th>(pretty)</th></tr></thead>' % extraColHeader)
        print '<tbody>'
        if len(times):
            lTimes = '%d' % len(times)
            mTimes = '%.2f' % mean(times)
            pmTimes = prettyTime(mean(times))
        else:
            lTimes = '0'
            mTimes = '--'
            pmTimes = '--'
        print('<tr style="background-color:#F6E8AE"><td>%s</td><td style="text-align:center;">%s</td>'
              '%s<td style="text-align:center;">%s</td><td style="text-align:center;">%s</td></tr>'
              % ('Overall', lTimes, extraCol, mTimes, pmTimes))
        for u in subList:
            if len(subTimesDict[u]):
                ulTimes = '%d' % len(subTimesDict[u])
                umTimes = '%.2f' % mean(subTimesDict[u])
                upmTimes = prettyTime(mean(subTimesDict[u]))
            else:
                ulTimes = '0'
                umTimes = '--'
                upmTimes = '--'
            print('<tr><td>%s</td><td style="text-align:center;">%s</td>'
                  '%s<td style="text-align:center;">%s</td><td style="text-align:center;">%s</td></tr>'
              % (u, ulTimes, extraCol, umTimes, upmTimes))
        
    else:
        if pre == 'Cycle':
            extraCol = '%20s' % '-- '
            extraColHeader = '%20s' % 'Mean Len. '
        else:
            extraCol = ' '
            extraColHeader = ' '
        print('%20s %4s %s%20s %20s' % ('Name', 'n', extraColHeader, 'Mean time (s)', '(pretty)'))
        print('%20s %4d %s%20.2f %20s'
              % ('Overall', len(times), extraCol, mean(times), prettyTime(mean(times))))
        for u in subList:
            print('%20s %4d %s%20.2f %20s'
              % (u, len(subTimesDict[u]), extraCol, mean(subTimesDict[u]), prettyTime(mean(subTimesDict[u]))))

    if pre != 'Cycle':
        if options.isHtml:
            print '</tbody></table>'
        return
    # print the chromosome times
    sortedChromTimes, chromTimesDict = getSortedChromTimesList(status)
    if options.isHtml:
        for c in sortedChromTimes:
            if len(c) > 22:
                chrom = ('<span title="%s" class="chrname">'
                         '%s ..(%d).. %s</span>' % (c, c[0:7], len(c) - 12, c[-7:]))
            else:
                chrom = '<span title="%s" class="chrname">%s</span>' % (c, c)
            if c in status.chromosomeLengthsDict:
                chromLengthStr = prettyLength(mean(getChrLenList(status.chromosomeLengthsDict, c)))
            else:
                chromLengthStr = ''
            print('<tr><td>%s</td><td style="text-align:center;">%d</td><td style="text-align:center;">%s</td>'
                  '<td style="text-align:center;">%.2f</td><td style="text-align:center;">%s</td></tr>'
                  % (chrom, len(chromTimesDict[c]),
                     chromLengthStr,
                     mean(chromTimesDict[c]), 
                     prettyTime(mean(chromTimesDict[c]))))
        print '</tbody></table>'
    else:
        for c in sortedChromTimes:
            if c in status.chromosomeLengthsDict:
                chromLengthStr = '%20s ' % prettyLength(mean(getChrLenList(status.chromosomeLengthsDict, c)))
            else:
                chromLengthStr = ' '
            print('%20s %4d %s%20.2f %20s'
              % (c, len(chromTimesDict[c]), chromLengthStr, 
                 mean(chromTimesDict[c]), prettyTime(mean(chromTimesDict[c]))))

def findStalledCycles(runDir, stepsDict, isHtml, htmlDir = ''):
    """findStalledCycles() calculates the average time cycles take and the
    standard deviation of that time and then it looks for cycles that are
    currently *in progress* and reports in progress cycles that are three
    SDs or more beyond the mean time.
    """
    n = 0
    for s in stepsDict:
        if stepsDict[s].complete:
            n += 1
    if n < 10:
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

def numRunningSteps(stepsDict):
    """ takes a stepsDict and returns the current 
    number of incomplete and running steps
    """
    n = 0
    for s in stepsDict:
        if not stepsDict[s].complete:
            if stepsDict[s].startTime != -1:
                n += 1
    return n

def listCurrentCycles(runDir, stepsDict, isHtml, options, htmlDir=''):
    """listCurretnCycles() lists all currently running cycles, their current
    step, current step runtime and current cycle runtime.
    """
    stepArray = [None]
    for c in ['Cycle', 'Stats']:
        for i in xrange(1,4):
            stepArray.append('%sStep%d_start' % (c, i))
            stepArray.append('%sStep%d_end' % (c, i))
    running = False
    if numRunningSteps(stepsDict):
        running = True
        if isHtml:
            print '<h3>Currently running cycles</h3>'
            print '<div style="margin-left:2em;">'
            print '<table cellpadding="5"><thead>'
            print ('<tr><th style="text-align:left;">%s</th><th>%s</th>'
                   '<th>%s</th><th>%s</th></tr></thead>' % ('Cycle', 'Total Time', 'Step', 'Step Time'))
            print '<tbody>'
        else:
            print 'Currently running cycles:'
            print '%30s %15s %25s %15s' % ('Cycle', 'Total Time', 'Step', 'Step Time')
    for p in stepsDict:
        if stepsDict[p].complete:
            continue
        if stepsDict[p].startTime == -1:
            continue
        runTime = time.time() - stepsDict[p].startTime
        stepName, stepTime = currentStepInfo(stepsDict[p], options.simDir)
        if isHtml:
            print('<tr><td>%s%s</a></td><td>%s</td><td>%25s</td>'
                  '<td>%s</td></tr>' 
                  % (str2link(p, htmlDir), p, prettyTime(runTime), stepName, prettyTime(stepTime)))
        else:
            print '%30s %15s %25s %15s' % (p, prettyTime(runTime), stepName, prettyTime(stepTime))
    if isHtml and running:
        print '</tbody></table></div>'

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
    if 'StatsendEpochUTC' in s.timeDict:
        stats = False
    if 'TransalignstartEpochUTC' in s.timeDict:
        transStep = 'TransalignStep'
        if stats:
            curStep += ' + ' + transStep
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
        print '<h3>Distribution of simulation step runtimes (%s)</h3>' % timeUnits
        print '<pre style="margin-left:2em;">'
    else:
        print 'Distribution of simulation step runtimes (%s):' % timeUnits
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
    returnItems.sort(reverse = True)
    if isHtml:
        print '<h3>List of step runtimes</h3>'
        print '<div style="margin-left:2em;">'
        print '<table cellpadding="5" border="1" bordercolor="#cccccc"><thead>'
        print '<tr><th>Time (s)</th><th>(pretty)</th><th>Cycle</th></tr></thead>'
        print '<tbody>'
    else:
        print 'List of step runtimes:'
        print '%12s %12s %30s' % ('Time (s)', '(Pretty)', 'Step name')
        print '%12s %12s %30s' % ('-' * 10, '-' * 10, '-' * 20)
    for k, v in returnItems:
        if isHtml:
            print('<tr><td style="text-align:right;">%10.2f</td>'
                  '<td style="text-align:center">%s</td><td>%s%s</a></td></tr>' 
                  % (float(k), prettyTime(k), str2link(v, htmlDir), v))
        else:
            print '%12.2f %12s %30s' % (float(k), prettyTime(k), v)
    if isHtml:
        print '</tbody></table></div>'
    
    
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
    if 'simNodeTree' not in status.variables:
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        if newickTree.iD == None:
            newickTree.iD = options.rootName
        status.simNodeTree = buildSimNodeTree(newickTree, options)
        status.variables['simNodeTree'] = True
        # printSimNodeTree(status.simNodeTree)

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
    if newickTree.iD == None:
            newickTree.iD = options.rootName
    status.numCompletedSteps, status.stepsDict = simStepUpdater(newickTree, options.stepLength,
                                                                status.stepsDict, 
                                                                status.simNodeTree, status, 
                                                                options)
    #####
    # now we find the longest continuous branch yet to be simulated
    # longBranchSteps is a Branch() object
    status.longBranchSteps = longestRemainingBranchFinder(status.simNodeTree, options.stepLength,
                                                          status.stepsDict, options) 
    if status.longBranchSteps.name == '':
        status.longBranchSteps.name = options.rootName
    if not status.longBranchSteps.longestPath:
        status.longBranchSteps.name = ''
    timeHandler(status, options)

def buildSimNodeTree(nt, options):
    if nt is None:
        return None
    root = SimNode()
    t = None
    while nt.distance > 0:
        nt.distance -= options.stepLength
        if nt.distance < 0:
            nt.distance = 0
        if t is None:
            root = SimNode()
            root.name = lsc.nameTree(nt)
            t = root
        else:
            p = t
            t = SimNode()
            p.left = t
            t.parent = p
            t.name = lsc.nameTree(nt)
    if t is None:
        root.name = lsc.nameTree(nt)
        t = root
    t.left = buildSimNodeTree(nt.left, options)
    if t.left is not None:
        t.left.parent = t
    t.right = buildSimNodeTree(nt.right, options)
    if t.right is not None:
        t.right.parent = t
    return root

def printSimNodeTree(sn):
    if sn is None:
        return
    if sn.parent is None:
        p = 'None'
    else:
        p = sn.parent.name
    if sn.left is None:
        l = 'None'
    else:
        l = sn.left.name
    if sn.right is None:
        r = 'None'
    else:
        r = sn.right.name
    print 'node: %s parent: %s left: %s right: %s' % (sn.name, p, l, r)
    printSimNodeTree(sn.left)
    printSimNodeTree(sn.right)

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
    f = open(filename, 'wb')
    cPickle.dump(status, f, 2) # 2 is the format protocol, 2 = binary
    f.close()
    
def printInfoTable(status, options):
    if options.isHtml:
        if 'isDone' not in status.variables:
            print '<form name="Reload">'
            print('<input type="checkbox" name="checkboxReload" onclick="StartTime()" checked="checked">%s'
                  % 'Autoreload (30s)')
            print '</form>'
        print '<h3>Information</h3>'
        print('<p style="margin-left:2em;">Generated at %s, %s</p>' 
              % (time.strftime("%a, %d %b %Y %H:%M:%S (%Z)", time.localtime()),
                 time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime())))
        print '<div style="margin-left:2em;">'
        print '<table cellpadding="5"><tr><td>\n'
        elmDiv = '</td><td>'
        rowDiv = '<td></tr><tr><td>'
    else:
        print 'Information'
        print 'Generated at %s, %s' % (time.strftime("%a, %d %b %Y %H:%M:%S (%Z)", time.localtime()),
                                       time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime()))
        elmDiv = ' '
        rowDiv = ' '
    
    info1 = ('Tot. tree len: %f (%d steps of %s)' 
             % (status.totalTreeLength, status.numTotalSteps,
                str(options.stepLength).rstrip('0')))
    if status.longBranchSteps.longestPath > 0.0:
        info2 = ('Longest remaining branch: %.1f %s' 
                 % (status.longBranchSteps.longestPath, status.workingCycleString))
    else:
        info2 = 'Longest remaining branch: --' 
    print '%s%s%s%s' % (info1,elmDiv, info2, rowDiv)
    info1 = ('Tot. stps taken: %d of %d (%2.2f%% complete)' 
             % (status.numCompletedSteps, status.numTotalSteps,
                100 * (float(status.numCompletedSteps) / float(status.numTotalSteps))))
    info2 = '%sElapsed CPU time: %s (ave: %s / step)' % (' ' * 5, status.elapsedTreeTimeStr, 
                                                         status.aveBranchTimeStr)
    print '%s%s%s%s' % (info1, elmDiv, info2, rowDiv)
                                                                                              
    if options.isHtml:
        status.remainingTimeStr = '<b>' + status.remainingTimeStr + '</b>'
    else:
        status.remainingTimeStr = status.remainingTimeStr
    info1 = 'ETtC: %31s ' % status.remainingTimeStr
    info2 = '%sEToC: %12s ' % (' ' * 4, status.estTimeOfCompStr)
    info3 = 'Elapsed wall-clock: %17s ' % prettyTime(status.elapsedTime)
    info4 = '%sETRL: %12s' % (' ' * 4, status.estTotalRunLength)
    print('%s%s%s%s\n%s%s%s\n' 
          % (info1, elmDiv, info2, rowDiv, info3, elmDiv, info4))
    if options.isHtml:
        print '</td></tr></table>'
        print '</div>'

def printTree(status, options):
    #####
    # Draw the Tree!
    if options.drawText:
        if options.isHtml:
            print '<h3>Simulation Tree Status</h3>'
            print '<pre style="margin-left:2em;">'
        nt = newickTreeParser(options.inputNewick, 0.0)
        drawText(nt, options, options.stepLength, status.totalTreeDepthSteps,
                  scale = options.scale, stepsDict = status.stepsDict,
                  isHtml = options.isHtml, directory = options.htmlDir)
    
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

def getSortedChromTimesList(status):
    # print the chromosome times
    chromTimesDict = {}
    for s in status.stepsDict:
        if 'chromosomes' not in status.stepsDict[s].timeDict:
            continue
        for c in status.stepsDict[s].timeDict['chromosomes']:
            if c not in chromTimesDict:
                chromTimesDict[c] = []
            chromTimesDict[c].append(status.stepsDict[s].timeDict['chromosomes'][c])
    return sorted(chromTimesDict, key = lambda c: mean(chromTimesDict[c]), reverse = True), chromTimesDict

def getChromTimesDictStep(status):
    chromTimesDictStep = {}
    for s in status.stepsDict:
        if 'chromosomes' not in status.stepsDict[s].timeDict:
            continue
        for c in status.stepsDict[s].timeDict['chromosomes']:
            if c not in chromTimesDictStep:
                chromTimesDictStep[c] = {}
            if not status.stepsDict[s].name in chromTimesDictStep[c]:
                chromTimesDictStep[c][status.stepsDict[s].name] = status.stepsDict[s].timeDict['chromosomes'][c]
            else:
                raise RuntimeError('sim step %s has identically named chromosomes "%s"' 
                                   % (name, c))
    return chromTimesDictStep

def getChrLenList(chrLenDict, c):
    l = []
    if c not in chrLenDict:
        return l
    for n in chrLenDict[c]:
        l.append(chrLenDict[c][n])
    return l

def printChrTimes(status, options):
    if not options.printChrTimes:
        return
    sortedChromTimes = getSortedChromTimesList(status)
    chromTimesDictStep = getChromTimesDictStep(status)
    if options.isHtml:
        print '<h3>Chromosome Lengths and Times</h3>'
        print '<table cellpadding="5" border="1" bordercolor="#cccccc"><thead>'
        print '<tr><th>Length (bp)</th><th>Time (s)</th></tr></thead>'
        print '<tbody>'
    else:
        print 'Chromosome Lengths and Times'
        print '#Length (bp)\tTime (s)'
    for c in sortedChromTimes:
        for n in status.chromosomeLengthsDict[c]:
            if options.isHtml:
                print ('<tr><td style="text-align:center;">%d</td>'
                       '<td style="text-align:center;">%d</td></tr>'
                       % (status.chromosomeLengthsDict[c][n], chromTimesDictStep[c][n]))
            else:
                print '%d\t%d' % (status.chromosomeLengthsDict[c][n], 
                                  chromTimesDictStep[c][n])
    if options.isHtml:
        print '</tbody></table>'

def main():
    usage = ('usage: %prog --simDir path/to/dir [options]\n\n'
             '%prog can be used to check on the status of a running or completed\n'
             'evolverSimControl simulation.')
    parser = OptionParser(usage = usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, parser)

    status = unpackData(os.path.join(options.simDir, '.simulationStatus.pickle'))
    
    if 'isDone' not in status.variables:
        collectData(options, status)

    if options.isHtml:
        pureHtmlStart()
        initHtml(status)
        
    printInfoTable(status, options)
    printTree(status, options)
    printStats(status, options)
    
    printChrTimes(status, options)
    
    if options.isHtml:
        finishHtml()

    if 'dontUpdate' not in status.variables:
        packData(status, os.path.join(options.simDir, '.simulationStatus.pickle'))

if __name__ == "__main__":
    main()
