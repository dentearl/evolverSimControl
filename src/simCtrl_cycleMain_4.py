#!/usr/bin/env python
"""
cycleMain_4.py
 dent earl, dearl (a) soe ucsc edu

cycleMain_4.py is a python wrapper for the evolver suite of genome
evolution tools. cycleMain_4 is the fourth in a series of four
wrappers that is written to interface with jobTree, a cluster
interface written by Benedict Paten. Other members of the wrappers
are:
    cycleMain_1.py cycleMain_2.py cycleMain_3.py cycleMain_4.py

cycleMain_4.py Handles:
 1) CDS align step.
 2)

7 October 2009
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
import xml.etree.ElementTree as ET
import os, subprocess, sys
from optparse import OptionParser
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimCycle   as LSY

programs = ['evolver_evo', 'evolver_cvt',
            'simCtrl_commandEval.py', 'evolver_gene_deactivate.sh',
            'simCtrl_cycleStats_1.py', 'simCtrl_completeTimestamp.py']
LSC.verifyPrograms(programs)
(EVO_BIN, CVT_BIN, CMD_EVAL_BIN, GDACT_BIN,
 STATS_BIN, TIMESTAMP_BIN) = programs

def usage():
    print "USAGE: %s --parent parentDir/ --child childDir --params globalParamsDir/ --jobFile JOB_FILE [optional: --step ]" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def main(argv):
    parser=OptionParser()
    LSC.standardOptions(parser)
    LSY.standardOptions(parser)
    (options, args) = parser.parse_args()
    LSC.standardOptionsCheck(options, usage)
    LSY.standardOptionsCheck(options, usage)
    LSC.subTypeTimestamp(os.path.join(options.childDir,'cycleInfo.xml'),
                         'main', 'cycleStep_3_cycleMain_3_end')
    LSC.subTypeTimestamp(os.path.join(options.childDir,'cycleInfo.xml'),
                         'main', 'cycleStep_4_cycleMain_4_start')

    xmlTree = ET.parse(options.jobFile)
    jobElm=xmlTree.getroot()
    childrenElm = xmlTree.find('children')

    ####################
    # gene decativation step
    gDActCMD  = CMD_EVAL_BIN
    gDActCMD += ' JOB_FILE "'
    CMD  = GDACT_BIN 
    CMD += ' ' + os.path.join( options.parentDir, 'annots.gff' )
    CMD += ' ' + os.path.join( options.childDir, 'intra', 'evannots.gff' )
    CMD += ' ' + os.path.join( options.childDir, 'annots.gff' )
    CMD += ' ' + EVO_BIN
    CMD += ' >& ' + os.path.join( options.childDir, 'logs', 'gene_deactivate.log' )
    gDActCMD += LSC.commandPacker( CMD ) + '"'
    
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = gDActCMD

    # Last step, add the cycle timestamp.
    followUpCommand  = CMD_EVAL_BIN
    followUpCommand += ' JOB_FILE "'
    followUpCommand += LSC.commandPacker( TIMESTAMP_BIN + ' --cycleDir ' + options.childDir + ' --timeType=main ' )
    CMD  = STATS_BIN
    CMD += ' --childDir ' + options.childDir
    CMD += ' --parentDir ' + options.parentDir
    if options.noMEs:
        CMD += ' --noMEs'
    CMD += ' --jobFile JOB_XML'
    followUpCommand += LSC.commandPacker( CMD )
    followUpCommand += '"'

    jobElm=xmlTree.getroot()
    jobElm.attrib['command'] = followUpCommand
    xmlTree.write( options.jobFile )

if __name__ == "__main__":
    main(sys.argv)
