#!/usr/bin/env python
"""
simTree.py
dent earl, dearl (a) soe ucsc edu
16 oct 2009
a recursive-ish script to control an entire tree's simulation.
Give the parent, and a newick tree as a command line
option and simTree launches an evolver cycle via jobTree.py.
"""
##############################
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import glob, os, sys, time
from datetime import datetime # for loggingo
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

programs = ['simCtrl_simTreeFollowUp.py', 'simCtrl_cycleMain_1.py',
            'simCtrl_commandEval.py']
LSC.verifyPrograms(programs)
(SIMTREE_FOLLOW_PY, CYCLEBEGIN_PY,
 CMD_EVAL_BIN) = programs

def usage():
    print 'USAGE: '+sys.argv[0]+' --parent <dir> --out [optional dir] --tree <newick tree in quotes> --params <parameter dir> --jobFile JOB_FILE '
    print __doc__
    sys.exit(2)

def main():
    parser=OptionParser()
    LSC.standardOptions(parser)
    LST.standardOptions(parser)
    (options, args) = parser.parse_args()
    LSC.standardOptionsCheck(options, usage)
    LST.standardOptionsCheck(options, usage)
    if(options.outDir != None):
        workingDir = options.outDir
    else:
        (workingDir,tail) = os.path.split(options.parentDir)
    if(len(glob.glob(options.parentDir+'/*rev')) == 0): 
        sys.stderr.write('%s: Error: no *.rev found in parent dir.\n' %(sys.argv[0]))
        usage()

    ##############################
    # End of pre-processing.
    ##############################
    newickTree = newickTreeParser(options.inputNewick, 0.0)
    xmlTree = ET.parse(options.jobFile)
    childrenElm = xmlTree.find('children')
    nextTree = newickTreeParser(options.inputNewick, 0.0)
    if newickTree.distance > 0:
        newChild = ET.SubElement(childrenElm, 'child')
        newChild.attrib['command'] = LST.cycleBranchCommandBuilder(newickTree, nextTree, 'stem', options.stepSize,
                                                                   workingDir, options.parentDir, options.gParamsDir,
                                                                   options.seed, options.logBranch, options.testTree)
    ##########
    # Follow up command
    followUpCommand = SIMTREE_FOLLOW_PY +\
                      ' --parent '+options.parentDir+\
                      ' --tree "'+options.inputNewick + '"'+\
                      ' --params '+options.gParamsDir +\
                      ' --step '+str(options.stepSize) +\
                      ' --seed '+options.seed +\
                      ' --jobFile JOB_FILE'
    if (options.outDir != None):
        followUpCommand = followUpCommand + ' --out ' + options.outDir
    if(options.saveParent):
        followUpCommand = followUpCommand + ' --saveParent '
    if(options.isContinue):
        followUpCommand = followUpCommand + ' --isContinue '
    if(options.isBranchChild):
        followUpCommand = followUpCommand + ' --isBranchChild '
    if(options.testTree):
        followUpCommand = followUpCommand + ' --testTree '
    if(options.logBranch):
        followUpCommand = followUpCommand + ' --logBranch '
        LST.branchLog( '%25s: %s\n' % ('simTreeFollowUp.py', followUpCommand))
    jobElm=xmlTree.getroot()
    jobElm.attrib['command'] = followUpCommand

    xmlTree.write(options.jobFile)


if __name__ == "__main__":
    main()
