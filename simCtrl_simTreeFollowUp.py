#!/usr/bin/env python
"""
 simTreeFollowUp.py
 dent earl, dearl (a) soe ucsc edu
 19 oct 2009
 a recursive-sh script to control an entire tree's simulation.
 Given the parent, and a newick tree as a command line
 option and simTreeFollowUp launches simTree.py as a jobTree.py
 child process (or two, if at a brach point).
"""
##############################
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import copy, os, shutil, sys, subprocess
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

programs = ['simCtrl_simTree.py']
LSC.verifyPrograms(programs)
SIMTREE_PY = programs[0]

def usage():
    sys.stderr.write('USAGE: '+sys.argv[0]+' --parent <dir> --tree <newick tree in quotes>'+\
          ' --params <parameter directory> --jobFile JOB_FILE ')
    sys.exit(2)
    
def main():
    parser=OptionParser()
    LSC.standardOptions(parser)
    LST.standardOptions(parser)
    (options, args) = parser.parse_args()
    LSC.standardOptionsCheck(options, usage)
    LST.standardOptionsCheck(options, usage)
    
    if(options.outDir != None):
        workingDir = os.path.abspath(options.outDir)
    else:
        (workingDir,tail) = os.path.split(options.parentDir)
    sys.stderr.write('followUp.py sees: %s' %(str(options.inputNewick)))
    newickTree = newickTreeParser(options.inputNewick, 0.0)
    if(newickTree.distance < options.stepSize):
        newickTree.distance = 0
    else:
        newickTree.distance = newickTree.distance - options.stepSize
    ##########
    # BRANCHES SHARE A COMMON PARENT, by Definition
    name=LST.nameTree(newickTree)
    commonParent = os.path.join(workingDir, name)
    
    xmlTree = ET.parse(options.jobFile)
    childrenElm = xmlTree.find('children')
    if(newickTree.distance <= 0):
        newickTree.distance = 0
        if(newickTree.internal):
            nextTree = newickTreeParser(options.inputNewick, 0.0)
            ####################
            # LEFT BRANCH
            newChild = ET.SubElement(childrenElm, 'child')
            newChild.attrib['command'] = LST.treeBranchCommandBuilder(newickTree.left, 'Left Branch', options,
                                                              commonParent, options.gParamsDir)
            ####################
            # RIGHT BRANCH
            newChild = ET.SubElement(childrenElm, 'child')
            newChild.attrib['command'] = LST.treeBranchCommandBuilder(newickTree.right, 'Right Branch', options,
                                                              commonParent, options.gParamsDir)
    else:
        ##########
        # STEM
        #####
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        if (newickTree.distance < options.stepSize):
            newickTree.distance = 0
        else:
            newickTree.distance = newickTree.distance - options.stepSize
        if(newickTree.distance > 0):
            newChild = ET.SubElement(childrenElm, 'child')
            newChild.attrib['command'] = LST.treeBranchCommandBuilder(newickTree, 'Stem', options, commonParent,
                                                              options.gParamsDir)

    newickTree = newickTreeParser(options.inputNewick, 0.0) # we check to see if this is a leaf.
    if((not options.saveParent) and (options.isContinue) and
       (not options.isBranchChild) and (newickTree.distance > 0)): # settings check
        shutil.rmtree(options.parentDir) # burn parent dir to the ground
        if(options.logBranch):
            (head, tail) = os.path.split(options.parentDir)
            LST.branchLog( '%25s: %s. saveParent:%d isContinue:%d distance:%f\n' % ('performing rmtree() on ',
                                                                                    tail, options.saveParent, options.isContinue, newickTree.distance))
    else:
        if(options.logBranch):
            (head, tail) = os.path.split(options.parentDir)
            LST.branchLog( 'will not delete parent (%s). saveParent:%d isContinue:%d isBranchChild:%d\n' % (tail, options.saveParent,
                                                                                                            options.isContinue, options.isBranchChild))
    jobElm=xmlTree.getroot()
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main()
