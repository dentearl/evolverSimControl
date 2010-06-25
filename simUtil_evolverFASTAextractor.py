#!/usr/bin/env python
"""
eval_evolverpairwiseMAFextractor.py
dent earl, dearl (a) soe ucsc edu
16 nov 2009
A script that will take the path to a simulation
out directory and makes repeated calls to cvt
(in fabulous PARALLEL-vision!) in order
to extract FASTA files from either all cycles, or
only the leaves.
"""
##############################
import glob
import os
import re
import subprocess
import sys
from optparse import OptionParser
from sonLib.bioio import newickTreeParser
import xml.etree.ElementTree as ET
import simulation.lib.libSimControl as LSC

programs = ['evolver_cvt', 'evolver_transalign', 'simUtil_fastaNameCorrector.py']
LSC.verifyPrograms(programs)
CVT_BIN     = programs[0]
TRANS_BIN   = programs[1]
NAMECORRECT = programs[2]

def usage():
    sys.stderr.write('USAGE: %s --simDir <dir> [optional: --leavesOnly --log --debug]\n' % (sys.argv[0]))
    sys.exit(2)

def logMessage(message):
    curr = os.curdir
    logPath = os.path.join(curr, 'pairwise_log.log')
    if(os.path.exists(logPath)):
        os.remove(logPath)
    FILE = open(logPath, 'w')
    FILE.write( '%s' % (message))
    FILE.close()

def initOptions(parser):
    parser.add_option('-i', '--simDir',dest='simDir',
                      help='Out directory from simTree.py.')
    parser.add_option('-l', '--logFile', dest='logFile',
                      help='Record output from evolver suite information.')
    parser.add_option('-s', '--leavesOnly', action='store_true', dest='leavesOnly',
                      default=False, help='Will exclude any directory with a "-" '
                      'in the name, which are taken to be internal nodes.')
    parser.add_option('-d', '--debug', action='store_true', dest='isDebug',
                      default=False, help='Won\'t issue jobs, just reports what would'
                      ' have happened')

def checkOptions(options):
    if (options.simDir == None):
        sys.stderr.write('%s: Error, specify --cycleDir.\n' % sys.argv[0])
        usage()
    if (not os.path.isdir(options.simDir)):
        sys.stderr.write('%s: Error, cycle directory "%s" not a directory!\n' % (sys.argv[0], options.simDir))
        usage()
    options.simDir = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        sys.stderr.write('%s: Error, unable to find simutaltionInfo.xml.\n' % sys.argv[0])
        usage()
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    treeObj = infoTree.find('rootDir')
    options.rootDir=treeObj.text

def directoriesOnly(list):
    """directoriesOnly() takes a list of items from a directory
    and purges anything from the list that is not a directory.
    """
    for i in list:
        if (not os.path.isdir(i)):
            list.remove(i)

def extractLeaves(nt, leafDict):
    """Given a newick tree object, it returns a dict of
    leaf objects. Operates recursively.
    """
    if nt == None:
        return None
    nt.distance=0
    if (nt.right == None) and (nt.left == None):
        leafDict[nt.iD] = 1
    else:
        extractLeaves(nt.right, leafDict=leafDict)
        extractLeaves(nt.left , leafDict=leafDict)    
    
def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    
    cycles = glob.glob(os.path.join(options.simDir, '*'))
    directoriesOnly(cycles)
    cmdArray=[]
    jobsArray=[]
    leaves = {}
    nt = newickTreeParser(options.inputNewick, 0.0)
    extractLeaves(nt, leaves)
    for i in cycles:
        if options.isDebug:
            print 'Now checking out cycle: %s' % os.path.basename(i)
        if(options.leavesOnly):
            if not os.path.basename(i) in leaves:
                if options.isDebug:
                    print 'Cycle is not a leaf: %s' % os.path.basename(i)
                continue
            if options.isDebug:
                print 'Cycle is a leaf: %s' % os.path.basename(i)
        nameA     = str(i.split('/')[-1:])
        nameA     = nameA.replace('[','')
        nameA     = nameA.replace(']','')
        cleanName = nameA.replace('\'','')
        cleanName = cleanName.replace('[','')
        cleanName = cleanName.replace(']','')
        inputRev  = os.path.join(i,'seq.rev')
        outputFa  = os.path.join(i,'seq.fa')
        outNameFa = os.path.join(i,'seq.name.fa')
        masked    = os.path.join(i,'seq.masked.fa')
        fastaCMD  = CVT_BIN +\
                    ' -fromrev '+inputRev+\
                    ' -tofasta '+outputFa+\
                    ' && '+\
                    NAMECORRECT +\
                    ' --name '+ cleanName+\
                    ' < '+outputFa+\
                    ' > '+outNameFa
        cmdArray.append(fastaCMD)
        if (options.logFile):
            logPath = os.path.join(os.curdir, options.logFile)
        else:
            logPath = os.path.abspath('/dev/null')
            logFILE = open(logPath, 'w')
        if options.isDebug:
            print 'Would have issued job: %s' % fastaCMD
            job = 'hi mom.'
        else:
            job = subprocess.Popen(fastaCMD, shell=True, stderr=logFILE, stdout=logFILE)
        jobsArray.append(job)
    #####
    # Make sure the jobs finished.
    if options.isDebug:
        return
    i=-1
    for job in jobsArray:
        job.wait()
        i=i+1
        if(job.returncode):
            if(options.log):
                logMessage('%s: Error in a sub-processes "%s", returncode: %s.\n' %(sys.argv[0], cmdArray[i], job.returncode))
            else:
                sys.stderr.write('%s: Error in a sub-processes "%s", returncode: %s.\n' %(sys.argv[0], cmdArray[i], job.returncode))
            

if __name__ == "__main__":
    main()
