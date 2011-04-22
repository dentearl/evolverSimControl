#!/usr/bin/env python
"""
cycleMain_1.py
 dent earl, dearl (a) soe ucsc edu

cycleMain_1.py is a python wrapper for the evolver suite of genome
evolution tools. cycleBegin is the first in a series of four
wrappers that is written to interface with jobTree, a cluster
interface written by Benedict Paten. Other members of the wrappers
are:
    cycleMain_1.py cycleMain_2.py cycleMain_3.py cycleMain_4.py

cycleMain_1.py Handles:
 1) Inter simulation step.

5 October 2009
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
import os, sys
import shutil # for rmtree, for the childDir
from optparse import OptionParser
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimCycle   as LSY

programs = ['evolver_evo', 'evolver_cvt',
            'simCtrl_cycleMain_2.py', 'simCtrl_commandEval.py',
            'evolver_handle_mobiles.pl']
LSC.verifyPrograms( programs )
( EVO_BIN, CVT_BIN, CYCLE_MAIN2, CMD_EVAL_BIN,
  MOBILES_BIN ) = programs

def usage():
    print "USAGE: %s --parent parentDir/ --child childDir --params globalParamsDir/ --jobFile JOB_FILE [optional: --step ]" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def buildDirs( options ):
    """
    bulidDirs accepts a string and creates the necessary
    directory structure for a cycle step.
    """
    try:
        os.mkdir(options.childDir)
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise
    else:
        # Step one, add the timestamp.
        os.mkdir(os.path.join(options.childDir,'intra'))
        os.mkdir(os.path.join(options.childDir,'inter'))
        os.mkdir(os.path.join(options.childDir,'logs'))
        if not options.noMEs:
            os.mkdir(os.path.join(options.childDir,'mobiles'))
        os.mkdir(os.path.join(options.childDir,'stats'))
        
        root=ET.Element('info')
        e=ET.SubElement(root, 'parentDir')
        e.text=options.parentDir
        e=ET.SubElement(root, 'childDir')
        e.text=options.childDir
        e=ET.SubElement(root, 'stepSize')
        e.text=str(options.stepSize).rstrip('0')
        e=ET.SubElement(root, 'cycleCommand')
        e.text = ' '.join(sys.argv)
        info=ET.ElementTree(root)
        info.write(os.path.join(options.childDir,'cycleInfo.xml'))
        LSC.typeTimestamp(os.path.join(options.childDir,'cycleInfo.xml'),
                          'main', 'Start')
        LSC.subTypeTimestamp(os.path.join(options.childDir,'cycleInfo.xml'),
                             'main', 'cycleStep_1_cycleMain_1_start')

def main(argv):
    parser = OptionParser()
    LSC.standardOptions( parser )
    LSY.standardOptions( parser )
    ( options, args ) = parser.parse_args()
    LSC.standardOptionsCheck( options, usage )
    LSY.standardOptionsCheck( options, usage )

    if os.path.exists( options.childDir ):
        sys.stderr.write( 'Directory [%s] already exists!\n' % options.childDir )
        sys.exit( 1 )

    # cycleBeginStep
    # first thing to do will be to make the directories
    buildDirs( options )
    
    # next will be the inter step
    interCMD = CMD_EVAL_BIN
    interCMD +=' JOB_FILE "'
    CMD  = EVO_BIN
    CMD += ' -interchr ' +os.path.join(options.parentDir, 'seq.rev')
    CMD += ' -inannots ' +os.path.join(options.parentDir, 'annots.gff')
    CMD += ' -aln ' + os.path.join(options.childDir, 'inter', 'inter.aln.rev')
    CMD += ' -outchrnames ' + os.path.join(options.childDir, 'inter', 'inter.chrnames.txt')
    CMD += ' -outannots ' + os.path.join(options.childDir, 'inter', 'inter.outannots.gff')
    CMD += ' -outseq ' + os.path.join(options.childDir, 'inter','inter.outseq.rev')
    CMD += ' -outgenome ' + options.theChild+'.inter'
    CMD += ' -branchlength ' + str(options.stepSize)
    CMD += ' -statsfile ' + os.path.join(options.childDir, 'stats', 'inter.stats')
    CMD += ' -model ' + os.path.join(options.gParamsDir, 'model.txt')
    CMD += ' -seed '+ options.seed
    CMD += ' -logevents '
    CMD += ' -log ' + os.path.join(options.childDir, 'logs', 'inter.log')
    ### CMD += ' -xfertandems' removed following may12 2010 crashes related to free() heap stack corruption
    interCMD +=LSC.commandPacker( CMD )
    interCMD += ' '
    if not options.noMEs:
        CMD  = MOBILES_BIN
        CMD += ' --evo '+EVO_BIN
        CMD += ' --cvt ' + CVT_BIN
        CMD += ' --py '+ os.path.dirname( EVO_BIN )
        CMD += ' --parentDir ' + options.parentDir
        CMD += ' --genome ' + options.theParent
        CMD += ' --stepSize '+ str( options.stepSize )
        CMD += ' --mefa ' + os.path.join( options.parentDir, 'mobiles', 'ME.fa' )
        CMD += ' --megff ' + os.path.join( options.parentDir, 'mobiles', 'ME.gff' )
        CMD += ' --ltr ' + os.path.join( options.parentDir, 'mobiles', 'LTR.fa' )
        CMD += ' --mescfg ' + os.path.join( options.gParamsDir, 'mes.cfg' )
        CMD += ' --model ' + os.path.join( options.gParamsDir, 'model.mes.txt' )
        CMD += ' > ' + os.path.join( options.childDir, 'logs', 'mobiles.log' )
        interCMD +=LSC.commandPacker( CMD )
        interCMD += ' '
        interCMD +=LSC.commandPacker('mv mes.fa '+os.path.join(options.childDir,'mobiles','mes.fa'))
        interCMD +=' '
        interCMD +=LSC.commandPacker('mv ME_output.fa '+os.path.join(options.childDir,'mobiles','ME.fa'))
        interCMD +=' '
        interCMD +=LSC.commandPacker('mv ME_output.gff '+os.path.join(options.childDir,'mobiles','ME.gff'))
        interCMD +=' '
        interCMD +=LSC.commandPacker('mv ME_output_ltrs.fa '+os.path.join(options.childDir,'mobiles','LTR.fa'))
    interCMD +='"'

    xmlTree = ET.parse(options.jobFile)
    jobElm=xmlTree.getroot()

    childrenElm = xmlTree.find('children')
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = interCMD
    
    # follow up job will be the intra steps (cycleMid_1.py)
    followUpCommand  = CYCLE_MAIN2
    followUpCommand += ' --parent ' + options.parentDir
    followUpCommand += ' --child '  + options.childDir
    followUpCommand += ' --step '   + str( options.stepSize )
    followUpCommand += ' --params ' + options.gParamsDir
    followUpCommand += ' --seed '   + options.seed
    if options.noMEs:
        followUpCommand += ' --noMEs'
    followUpCommand += ' --jobFile JOB_FILE '
    jobElm.attrib['command'] = followUpCommand
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)
