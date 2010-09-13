#!/usr/bin/env python
"""
postSimRepeatMaskWrapper.py
dent earl, dearl (a) soe ucsc edu
9 september 2010
a script to be run following the completion of a
simulation. The script will check the simulation dir
and then run the repeatMasking_doCluster.py on each
of the cycle/seq.name.fa files.
"""
##############################
import glob
import os
import subprocess
import sys
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

programs = ['repeatMasking_doCluster.py', 'parasol']
LSC.verifyPrograms(programs)
( REPEATMASK_BIN, PARASOL ) = programs

def usage():
    print 'USAGE: '+sys.argv[0]+' --simDir <dir> [ optional: --maxJobs <number> --log ]'
    print __doc__
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-p', '--simDir',dest='simDir',
                      help='Simulation directory.')
    parser.add_option('-l', '--log', dest='logFile',
                      default=False, action='store_true',
                      help='Record output from evolver suite information.')
    parser.add_option('-w', '--wait', dest='isWait',
                      default=False, action='store_true',
                      help='Holds the script open until all the cluster jobs return. For use in a batch script.')
    parser.add_option('-m', '--maxJobs',dest='maxJobs',
                      type='int', default=200,
                      help='Total number of max jobs to divide among all repeat masking instances.')

def checkOptions(options):
    if (options.simDir == None):
        sys.stderr.write('%s: Error, specify simulation dir.\n' % sys.argv[0])
        usage()
    options.simDir  = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        sys.stderr.write('%s: Error, unable to find simulationInfo.xml.\n' % sys.argv[0])
        usage()
    if options.maxJobs < 1:
        sys.stderr.write('%s: Error, --maxJobs of %d too small.\n' % (sys.argv[0], options.maxJobs))
        usage()
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    treeObj = infoTree.find('rootDir')
    options.rootDir=treeObj.text

def findDirectories( simDir ):
    allItems  = glob.glob(os.path.join(simDir, '*'))
    allDirs   = []
    readyDirs = []
    for a in allItems:
        if os.path.isdir(a):
            allDirs.append(a)
    for a in allDirs:
        if os.path.exists(os.path.join(a,'seq.name.fa')):
            readyDirs.append(a)
    return readyDirs

def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    dirsToProcess = findDirectories( options.simDir )
    jobsArray = [] # not used, here for future use
    cmdArray  = []
    numJobs = min(1, int( options.maxJobs / len( dirsToProcess ) ))
    for thisDir in dirsToProcess:
        if (options.logFile):
            logPath = os.path.abspath(os.path.join(thisDir, 'logs', 'repeatMasking.log'))
        else:
            logPath = os.path.abspath('/dev/null')
        logFILE = open(logPath, 'w')
        maskCMD =  REPEATMASK_BIN
        maskCMD += ' --genome '  + os.path.join(thisDir, 'seq.name.fa')
        maskCMD += ' --workDir ' + os.path.join(thisDir, 'repeatMask')
        maskCMD += ' --maxJob=' + str(numJobs)
        cmdArray.append(maskCMD)
        job = subprocess.Popen(maskCMD, shell=True, stderr=logFILE, stdout=logFILE)
        jobsArray.append(job)
    if options.isWait:
        i=-1
        for job in jobsArray:
            i=i+1
            job.wait()
            if(job.returncode):
                sys.stderr.write('%s: Error in a sub-processes "%s", returncode: %s.\n' %(sys.argv[0], cmdArray[i], job.returncode))
        
if __name__ == "__main__":
    main()
