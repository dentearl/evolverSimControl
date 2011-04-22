#!/usr/bin/env python
"""
simCtrl_checkCycleTimes.py
dent earl, dearl (a) soe ucsc edu
29 June 2010

This script will look at a simCtrl output directory
and will report back all the runtimes for the completed
cycles along with their tree-depth level.

I believe there may be some correlation between
tree-depth and runtime and I want to verify this.

"""
##############################
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

def usage():
    print 'USAGE: %s --help' %(sys.argv[0])
    print __doc__
    sys.exit(2)

class Cycle:
    """I use 'Cycle' objects to keep track of cycle
    names, tree-depth, and runtimes.
    """
    def __init__(self):
        self.name=''
        self.nameWithoutDistance=''
        self.depth=0
        self.time=0.0

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

def initOptions(parser):
    parser.add_option('-d', '--dir',dest='dir',
                      help='Parent directory.')
    parser.add_option('--formatR',dest='isR', action='store_true',default=False,
                      help='prints output in an R friendly way')

def checkOptions(options):
    options.dir=os.path.abspath(options.dir)
    if not os.path.exists(os.path.join(options.dir, 'simulationInfo.xml')):
        sys.stderr.write('%s: Error, unable to locate simulationInfo.xml.\n' % sys.argv[0])
        usage()
    infoTree = ET.parse(os.path.join(options.dir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick = treeObj.text
    nt= newickTreeParser(options.inputNewick, 0.0)
    if nt.distance==0:
        if nt.iD == None:
            nt.iD = 'root'
            options.inputNewick = printBinaryTree(nt,1)
            
    treeObj = infoTree.find('stepSize')
    stepSize=treeObj.text
    options.stepSize=float(stepSize)

def directoriesOnly(list):
    """directoriesOnly() takes a list of items from a directory
    and purges anything from the list that is not a directory.
    """
    for i in list:
        if (not os.path.isdir(i)):
            list.remove(i)
            
def extractCycles(options, status):
    """extract cycles should return a list of Cycle objects for
    the simulation.
    """
    cycleList = []
    # just an fyi, in python .append on lists *is* atomic.
    nt = newickTreeParser(options.inputNewick, 1.0)
    recursiveCycleExtraction(nt, options, cycleList)
    status.cycleList = cycleList

def recursiveCycleExtraction(nt, options, cycleList, depth=1):
    if (nt == None):
        return
    nt.distance = nt.distance - options.stepSize
    if nt.distance < 0:
        nt.distance = 0
    name = LST.nameTree(nt)
    nameWOD = LST.nameTree(nt, reportDistance=0)
    c = Cycle()
    c.name = name
    c.nameWithoutDistance = nameWOD
    c.depth = depth
    t = getCycleRunTime(os.path.join(options.dir, name, 'cycleInfo.xml'))
    if t != 0.0:
        c.time = t
        cycleList.append(c)
    if nt.distance == 0:
        recursiveCycleExtraction(nt.right, options, cycleList, depth + 1 )
        recursiveCycleExtraction(nt.left, options, cycleList, depth + 1 )
    else:
        recursiveCycleExtraction(nt, options, cycleList, depth + 1 )
    
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

def printCycles(options, status):
    if options.isR:
        sys.stdout.write('depth = c(')
        for i in range(len(status.cycleList)):
            if i == len(status.cycleList)-1:
                s = ')\n'
            else:
                s = ', '
            sys.stdout.write( "%2d%s" %(status.cycleList[i].depth, s))
        sys.stdout.write('time = c(')
        for i in range(len(status.cycleList)):
            if i == len(status.cycleList)-1:
                s = ')\n'
            else:
                s = ', '
            sys.stdout.write( "%f%s" %(status.cycleList[i].time, s))
        sys.stdout.write('name = c(')
        for i in range(len(status.cycleList)):
            if i == len(status.cycleList)-1:
                s = ')\n'
            else:
                s = ', '
            sys.stdout.write( "'%s'%s" %(status.cycleList[i].nameWithoutDistance, s))
        print 'df = data.frame(Time=time, Depth=depth, Name=name)'
        print 'require(ggplot2)'
        print 'qplot(x=Depth, y=Time, data=df, color=Name, xlab="Tree Depth", ylab="Time (sec)")'
    else:
        for c in status.cycleList:
            print "%2d %10.2f %s" %(c.depth, c.time, c.name)

def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    status = Status()
    status.cycleDirs = glob.glob(os.path.join(options.dir,'*'))
    directoriesOnly(status.cycleDirs)
    extractCycles(options, status)
    printCycles(options, status)

if __name__ == "__main__":
    main()
