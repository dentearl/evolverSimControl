#!/usr/bin/env python
"""
cycleMain_1.py
 dent earl, dearl (a) soe ucsc edu

cycleMain_1.py is a python wrapper for the evolver suite of genome
evolution tools. cycleBegin is the first in a series of four
wrappers that is written to interface with jobTree.py, a cluster
interface written by Benedict Paten. Other members of the wrappers
are:
    cycleMain_1.py cycleMain_2.py cycleMain_3.py cycleMain_4.py

cycleMain_1.py Handles:
 1) Inter simulation step.

5 October 2009
"""
######################################
import xml.etree.ElementTree as ET
import os, sys
import shutil # for rmtree, for the childDir
from optparse import OptionParser
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimCycle   as LSY

programs = ['evolver_evo', 'evolver_cvt',
            'simCtrl_cycleMain_2.py', 'simCtrl_commandEval.py',
            'evolver_handle_mobiles.pl']
LSC.verifyPrograms(programs)
(EVO_BIN, CVT_BIN, CYCLE_MID1, CMD_EVAL_BIN,
 MOBILES_BIN) = programs

def usage():
    print "USAGE: %s --parent parentDir/ --child childDir --params globalParamsDir/ --jobFile JOB_FILE [optional: --step ]" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def buildDirs(options):
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
    parser=OptionParser()
    LSC.standardOptions(parser)
    LSY.standardOptions(parser)
    (options, args) = parser.parse_args()
    LSC.standardOptionsCheck(options, usage)
    LSY.standardOptionsCheck(options, usage)

    if os.path.exists(options.childDir):
        shutil.rmtree(options.childDir)   # yup, burn it to the ground! That'll learn ya to have name collisions!

    # cycleBeginStep
    # first thing to do will be to make the directories
    buildDirs(options)
    
    # next will be the inter step
    interCMD = CMD_EVAL_BIN
    interCMD +=' JOB_FILE "'
    interCMD +=LSC.commandPacker(EVO_BIN +\
                                 ' -interchr ' +os.path.join(options.parentDir, 'seq.rev') + \
                                 ' -inannots ' +os.path.join(options.parentDir, 'annots.gff') + \
                                 ' -aln ' + os.path.join(options.childDir, 'inter', 'inter.aln.rev') + \
                                 ' -outchrnames ' +  os.path.join(options.childDir, 'inter', 'inter.chrnames.txt') + \
                                 ' -outannots ' +  os.path.join(options.childDir, 'inter', 'inter.outannots.gff') + \
                                 ' -outseq ' + os.path.join(options.childDir, 'inter','inter.outseq.rev') + \
                                 ' -outgenome ' +  options.theChild+'.inter' + \
                                 ' -branchlength '  +str(options.stepSize) + \
                                 ' -statsfile ' + os.path.join(options.childDir, 'stats', 'inter.stats') + \
                                 ' -model ' + os.path.join(options.gParamsDir, 'model.txt') +\
                                 ' -seed '+ options.seed+\
                                 ' -logevents ' + \
                                 ' -log ' + os.path.join(options.childDir, 'logs', 'inter.log')) # ' -xfertandems' removed following may12 crashes related to free() heap stack corruption
    interCMD += ' '
    interCMD +=LSC.commandPacker(MOBILES_BIN+\
                                 ' --evo '+EVO_BIN+\
                                 ' --cvt '+CVT_BIN+\
                                 ' --py '+ os.path.dirname(EVO_BIN)+\
                                 ' --parentDir '+options.parentDir+\
                                 ' --genome '+options.theParent+\
                                 ' --stepSize '+str(options.stepSize)+\
                                 ' --mefa '+os.path.join(options.parentDir,'mobiles', 'ME.fa')+\
                                 ' --megff '+os.path.join(options.parentDir,'mobiles', 'ME.gff')+\
                                 ' --ltr '+os.path.join(options.parentDir,'mobiles', 'LTR.fa')+\
                                 ' --mescfg '+os.path.join(options.gParamsDir,'mes.cfg')+\
                                 ' --model '+os.path.join(options.gParamsDir,'model.mes.txt')+\
                                 ' > '+os.path.join(options.childDir,'logs','mobiles.log'))
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
    followUpCommand = CYCLE_MID1 + ' --parent ' + options.parentDir +\
                      ' --child '  + options.childDir +\
                      ' --step '   + str(options.stepSize) +\
                      ' --params ' + options.gParamsDir +\
                      ' --seed '   + options.seed+\
                      ' --jobFile JOB_FILE '
    jobElm.attrib['command'] = followUpCommand
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)
