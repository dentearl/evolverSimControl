#!/usr/bin/env python
"""
cycleMain_2.py, 19 October 2009
 dent earl, dearl (a) soe ucsc edu

cycleMain_2.py is a python wrapper for the evolver suite of genome
evolution tools. cycleMain_2 is the second in a series of four 
wrappers that is written to interface with jobTree, a cluster
interface written by Benedict Paten. Other members of the wrappers
are:
    cycleMain_1.py cycleMain_2.py cycleMain_3.py cycleMain_4.py 

cycleMain_2.py Handles:
 1) 1st Transalign step.
 2) All intra simulation.
*3) Tandem Repeats Finder

* New features Feb 2010
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
import sys, os
from optparse import OptionParser
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimCycle   as LSY

programs = ['evolver_evo', 'evolver_cvt', 'simCtrl_cycleMain_3.py',
            'simCtrl_wrapperTRF.py', 'mv' , 'echo', 'sleep',
            'simCtrl_commandEval.py']
LSC.verifyPrograms(programs)
(EVO_BIN, CVT_BIN, CYCLE_MAIN3, TRF_WRAP_BIN,
 MV_BIN, ECHO_BIN, SLEEP_BIN, CMD_EVAL_BIN) = programs

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
                         'main', 'cycleStep_1_cycleMain_1_end')
    LSC.subTypeTimestamp(os.path.join(options.childDir,'cycleInfo.xml'),
                         'main', 'cycleStep_2_cycleMain_2_start')
    xmlTree = ET.parse(options.jobFile)
    jobElm=xmlTree.getroot()
    childrenElm = xmlTree.find('children')
    
    ##############################
    # Intra steps:
    # all chromosomes should be launched in child processes
    FILE = open(os.path.join(options.childDir, 'inter', 'inter.chrnames.txt'))
    for line in FILE:
        chrom = line.rstrip()
        newChild = ET.SubElement(childrenElm, 'child')
        intraCMD  = CMD_EVAL_BIN
        intraCMD += ' --statXML '+os.path.join(options.childDir, 'logs','micro.'+chrom+'.info.xml')
        intraCMD += ' JOB_FILE "'
        CMD  = EVO_BIN
        CMD += ' -inseq '    +os.path.join(options.childDir, 'inter', 'inter.outseq.rev')
        CMD += ' -chrname '  +chrom
        CMD += ' -branchlength '  +str(options.stepSize)
        CMD += ' -seed '     +options.seed
        if not options.noMEs:
            CMD += ' -mes '      +os.path.join(options.childDir, 'mobiles', 'mes.fa')
        CMD += ' -inannots ' +os.path.join(options.childDir, 'inter', 'inter.outannots.gff')
        CMD += ' -statsfile '+os.path.join(options.childDir, 'stats', chrom+'.stats')
        CMD += ' -codonsubs '+os.path.join(options.childDir, 'intra', chrom+'.codonsubs')
        CMD += ' -outannots '+os.path.join(options.childDir, 'intra', chrom+'.outannots.gff')
        CMD += ' -outgenome '+options.theChild
        CMD += ' -model '    +os.path.join(options.gParamsDir,'model.txt')
        CMD += ' -aln '      +'LOCAL_DIR/'+chrom+'.aln.rev'
        CMD += ' -outseq '   +'LOCAL_DIR/'+chrom+'.outseq.rev'
        CMD += ' -log '      +os.path.join(options.childDir,  'logs', 'intra.'+chrom+'.log')
        intraCMD += LSC.commandPacker( CMD )
        
        CMD  = CVT_BIN
        CMD += ' -fromrev ' + 'LOCAL_DIR/' + chrom + '.outseq.rev'
        CMD += ' -tofasta ' + 'LOCAL_DIR/' + chrom + '.outseq.fa'
        CMD += ' -log '     + os.path.join( options.childDir,  'logs', 'intra.'+chrom+'.tofasta.log' )
        intraCMD += LSC.commandPacker( CMD )
        
        CMD  = TRF_WRAP_BIN + ' LOCAL_DIR/'+chrom+'.outseq.fa 2 7 7 80 10 50 500 -d -h '
        intraCMD += LSC.commandPacker( CMD )
        intraCMD += LSC.commandPacker( ECHO_BIN + ' \' \' ' )
        intraCMD += LSC.commandPacker( MV_BIN + ' LOCAL_DIR/' + chrom + '.*.dat ' + os.path.join( options.childDir, 'intra/' ) )
        intraCMD += LSC.commandPacker( MV_BIN + ' LOCAL_DIR/' + chrom + '.outseq.rev ' + os.path.join( options.childDir, 'intra/' ) + ' ' )
        intraCMD += LSC.commandPacker( MV_BIN + ' LOCAL_DIR/' + chrom + '.aln.rev '+ os.path.join( options.childDir, 'intra/' ) + ' ' )
        intraCMD +='"'
        newChild.attrib['command'] = intraCMD
    FILE.close()
    
    jobElm=xmlTree.getroot()
    followUpCommand  = CYCLE_MAIN3
    followUpCommand += ' --parent ' + options.parentDir
    followUpCommand += ' --child '  + options.childDir
    followUpCommand += ' --params ' + options.gParamsDir
    followUpCommand += ' --seed '   + options.seed
    if options.noMEs:
        followUpCommand += ' --noMEs'
    followUpCommand += ' --jobFile JOB_FILE '
    jobElm.attrib['command'] = followUpCommand
    ##############################
    # finished
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)
