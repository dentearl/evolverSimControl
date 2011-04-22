#!/usr/bin/env python
"""
simTree.py
dent earl, dearl (a) soe ucsc edu
16 oct 2009
a recursive-ish script to control an entire tree's simulation.
Give the parent, and a newick tree as a command line
option and simTree launches an evolver cycle via jobTree.
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
import glob, os, sys, time
from datetime import datetime # for loggingo
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

programs = [ 'simCtrl_simTreeFollowUp.py', 'simCtrl_cycleMain_1.py',
             'simCtrl_commandEval.py' ]
LSC.verifyPrograms( programs )
( SIMTREE_FOLLOW_PY, CYCLEBEGIN_PY,
  CMD_EVAL_BIN ) = programs

def usage():
    print 'USAGE: '+sys.argv[0]+' --parent <dir> --out [optional dir] --tree <newick tree in quotes> --params <parameter dir> --jobFile JOB_FILE '
    print __doc__
    sys.exit(2)

def main():
    parser=OptionParser()
    LSC.standardOptions( parser )
    LST.standardOptions( parser )
    ( options, args ) = parser.parse_args()
    LSC.standardOptionsCheck( options, usage )
    LST.standardOptionsCheck( options, usage )
    if options.outDir != None:
        workingDir = options.outDir
    else:
        ( workingDir, tail ) = os.path.split( options.parentDir )
    if len( glob.glob( options.parentDir+'/*rev' ) ) == 0: 
        sys.stderr.write('%s: Error: no *.rev found in parent dir.\n' %( sys.argv[0] ))
        usage()

    ##############################
    # End of pre-processing.
    ##############################
    newickTree = newickTreeParser( options.inputNewick, 0.0 )
    xmlTree = ET.parse(options.jobFile)
    childrenElm = xmlTree.find('children')
    nextTree = newickTreeParser( options.inputNewick, 0.0 )
    if newickTree.distance > 0:
        newChild = ET.SubElement( childrenElm, 'child' )
        newChild.attrib['command'] = LST.cycleBranchCommandBuilder( newickTree, nextTree, 'stem', options.stepSize,
                                                                    workingDir, options.parentDir, options.gParamsDir,
                                                                    options.seed, options.logBranch, options.testTree,
                                                                    options.noMEs )
    ##########
    # Follow up command
    followUpCommand  = SIMTREE_FOLLOW_PY
    followUpCommand += ' --parent ' + options.parentDir
    followUpCommand += ' --tree "' + options.inputNewick + '"'
    followUpCommand += ' --params ' + options.gParamsDir
    followUpCommand += ' --step ' + str( options.stepSize )
    followUpCommand += ' --seed '+ options.seed 
    followUpCommand += ' --jobFile JOB_FILE'
    if options.outDir != None:
        followUpCommand = followUpCommand + ' --out ' + options.outDir
    if options.removeParent:
        followUpCommand = followUpCommand + ' --removeParent'
    if options.isContinue:
        followUpCommand = followUpCommand + ' --isContinue'
    if options.isBranchChild:
        followUpCommand = followUpCommand + ' --isBranchChild'
    if options.testTree:
        followUpCommand = followUpCommand + ' --testTree'
    if options.noMEs:
        followUpCommand = followUpCommand + ' --noMEs'
    if options.logBranch:
        followUpCommand = followUpCommand + ' --logBranch'
        LST.branchLog( '%25s: %s\n' % ('simTreeFollowUp.py', followUpCommand))
    jobElm=xmlTree.getroot()
    jobElm.attrib['command'] = followUpCommand

    xmlTree.write(options.jobFile)


if __name__ == "__main__":
    main()
