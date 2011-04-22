#!/usr/bin/env python
"""
 simTreeFollowUp.py
 dent earl, dearl (a) soe ucsc edu
 19 oct 2009
 a recursive-sh script to control an entire tree's simulation.
 Given the parent, and a newick tree as a command line
 option and simTreeFollowUp launches simTree.py as a jobTree
 child process (or two, if at a brach point).
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
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
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
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import copy, os, shutil, sys, subprocess
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

programs = ['simCtrl_simTree.py', 'simCtrl_cycleSerialTransalign.py']
LSC.verifyPrograms(programs)
(SIMTREE_PY, CYCLETRANS_PY) = programs

def usage():
    sys.stderr.write('USAGE: '+sys.argv[0]+' --parent <dir> --tree <newick tree in quotes>'+\
          ' --params <parameter directory> --jobFile JOB_FILE ')
    sys.exit(2)

def cycleIsLeaf(nt):
    """Returns True if the newick tree supplied has 0 distance
    and is not a node
    """
    if (not nt.internal) and (nt.distance == 0):
        return True
    return False
    
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
        if newickTree.internal:
            nextTree = newickTreeParser(options.inputNewick, 0.0)
            ####################
            # LEFT BRANCH
            if not cycleIsLeaf(newickTree.left):
                newChild = ET.SubElement(childrenElm, 'child')
                newChild.attrib['command'] = LST.treeBranchCommandBuilder(newickTree.left, 'Left Branch', options,
                                                                          commonParent, options.gParamsDir)
                # Left branches look backwards for the transalign,
                # right braches do not.
                name = LST.nameTree(newickTree.left)
                childDir = os.path.join(workingDir, name)
                newChild = ET.SubElement(childrenElm, 'child')
                transCMD = CYCLETRANS_PY +\
                           ' --targetDir ' + commonParent +\
                           ' --jobFile JOB_FILE'
                newChild.attrib['command'] = transCMD
                LST.commandRecorder(transCMD, commonParent)
            else:
                # cycle is a leaf
                name = LST.nameTree(newickTree.left)
                childDir = os.path.join(workingDir, name)
                followUpCommand = CYCLETRANS_PY +\
                                  ' --targetDir ' + commonParent +\
                                  ' --jobFile JOB_FILE'
                jobElm=xmlTree.getroot()
                jobElm.attrib['command'] = followUpCommand
                LST.commandRecorder(followUpCommand, childDir)
            ####################
            # RIGHT BRANCH
            if not cycleIsLeaf(newickTree.left):
                newChild = ET.SubElement(childrenElm, 'child')
                newChild.attrib['command'] = LST.treeBranchCommandBuilder(newickTree.right, 'Right Branch', options,
                                                                          commonParent, options.gParamsDir)
            else:
                # cycle is a leaf
                name = LST.nameTree(newickTree.right)
                childDir = os.path.join(workingDir, name)
                followUpCommand = CYCLETRANS_PY +\
                                  ' --targetDir ' + commonParent +\
                                  ' --jobFile JOB_FILE'
                jobElm=xmlTree.getroot()
                jobElm.attrib['command'] = followUpCommand
                LST.commandRecorder(followUpCommand, childDir)
        else:
            ##########
            # STEM without distance that isn't internal (i.e., a leaf!)
            if cycleIsLeaf(newickTree):
                name = LST.nameTree(newickTree)
                childDir = os.path.join(workingDir, name)
                followUpCommand = CYCLETRANS_PY +\
                                  ' --targetDir ' + commonParent +\
                                  ' --jobFile JOB_FILE'
                jobElm=xmlTree.getroot()
                jobElm.attrib['command'] = followUpCommand
                LST.commandRecorder(followUpCommand, childDir)
    else:
        ##########
        # STEM with distance
        #####
        newickTree = newickTreeParser(options.inputNewick, 0.0)
        if newickTree.distance < options.stepSize:
            newickTree.distance = 0
        else:
            newickTree.distance = newickTree.distance - options.stepSize
        if newickTree.distance > 0:
            newChild = ET.SubElement(childrenElm, 'child')
            newChild.attrib['command'] = LST.treeBranchCommandBuilder(newickTree, 'Stem', options, commonParent,
                                                                      options.gParamsDir)
            name = LST.nameTree(newickTree.left)
            childDir = os.path.join(workingDir, name)
            newChild = ET.SubElement(childrenElm, 'child')
            transCMD = CYCLETRANS_PY +\
                       ' --targetDir ' + commonParent +\
                       ' --jobFile JOB_FILE'
            newChild.attrib['command'] = transCMD
            LST.commandRecorder(transCMD, commonParent)

    newickTree = newickTreeParser(options.inputNewick, 0.0) # we check to see if this is a leaf.
    if( (options.removeParent) and (options.isContinue) and
       (not options.isBranchChild) and (newickTree.distance > 0) ): # settings check
        shutil.rmtree(options.parentDir) # burn parent dir to the ground
        if options.logBranch:
            (head, tail) = os.path.split(options.parentDir)
            LST.branchLog( '%25s: %s. removeParent:%d isContinue:%d distance:%f\n' % ('performing rmtree() on ',
                                                                                      tail, options.removeParent, options.isContinue, newickTree.distance))
    else:
        if options.logBranch:
            (head, tail) = os.path.split(options.parentDir)
            LST.branchLog( 'will not delete parent (%s). removeParent:%d isContinue:%d isBranchChild:%d\n' % (tail, options.removeParent,
                                                                                                              options.isContinue, options.isBranchChild))
    jobElm=xmlTree.getroot()
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main()
