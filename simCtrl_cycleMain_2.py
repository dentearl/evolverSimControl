#!/usr/bin/env python
"""
cycleMain_2.py, 19 October 2009
 dent earl, dearl (a) soe ucsc edu

cycleMain_2.py is a python wrapper for the evolver suite of genome
evolution tools. cycleMain_2 is the second in a series of four 
wrappers that is written to interface with jobTree.py, a cluster
interface written by Benedict Paten. Other members of the wrappers
are:
    cycleMain_1.py cycleMain_2.py cycleMain_3.py cycleMain_4.py 

cycleMain_2.py Handles:
 1) 1st Transalign step.
 2) All intra simulation.
*3) Tandem Repeats Finder

* New features Feb 2010
"""
########################################
import xml.etree.ElementTree as ET
import sys, os
from optparse import OptionParser
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimCycle   as LSY

programs = ['evolver_evo', 'evolver_cvt', 'evolver_transalign',
            'simCtrl_cycleMain_3.py', 'ln', 'simCtrl_wrapperTRF.py', 
            'mv' , 'echo', 'sleep', 'simCtrl_commandEval.py']
LSC.verifyPrograms(programs)
(EVO_BIN, CVT_BIN, TRANS_BIN, CYCLE_MAIN3, LINK_BIN, TRF_WRAP_BIN,
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
    
    # 1st transalign
    if( os.path.isfile(os.path.join(options.parentDir, 'root.aln.rev'))):
        # align the inter step to the root via the parent's 'to-root' alignment.
        transCMD = CMD_EVAL_BIN+\
                   ' --statXML '+os.path.join(options.childDir, 'logs','trans.1.info.xml')+\
                   ' JOB_FILE "'
        transCMD +=LSC.commandPacker(TRANS_BIN+\
                                     ' -in1 '+ os.path.join(options.childDir, 'inter','inter.aln.rev') + \
                                     ' -in2 '+ os.path.join(options.parentDir,'root.aln.rev')  + \
                                     ' -out '+ 'LOCAL_DIR/'+options.theParent+'.inter.aln.rev '+\
                                     ' -log '+ 'LOCAL_DIR/transalign.inter.log')
        transCMD +=LSC.commandPacker(MV_BIN+\
                                     ' '+'LOCAL_DIR/'+options.theParent+'.inter.aln.rev '+ os.path.join(options.childDir, 'inter/'))
        transCMD +=LSC.commandPacker(MV_BIN+\
                                     ' '+'LOCAL_DIR/transalign.inter.log '+ os.path.join(options.childDir, 'logs/'))
        transCMD +='"'
    else:
        # base case, the parent *is* the root.
        transCMD = CMD_EVAL_BIN+\
                   ' JOB_FILE "'+\
                   LSC.commandPacker(LINK_BIN +\
                                     ' -s inter.aln.rev '+os.path.join(options.childDir,'inter', options.theParent+'.inter.aln.rev'))+\
                   '"'

    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = transCMD

    ##############################
    # Intra steps:
    # all chromosomes should be launched in child processes
    FILE = open(os.path.join(options.childDir, 'inter', 'inter.chrnames.txt'))
    for line in FILE:
        chrom = line.rstrip()
        newChild = ET.SubElement(childrenElm, 'child')
        intraCMD = CMD_EVAL_BIN+\
                   ' --statXML '+os.path.join(options.childDir, 'logs','micro.'+chrom+'.info.xml')+\
                   ' JOB_FILE "'
        intraCMD +=LSC.commandPacker(EVO_BIN +\
                                     ' -inseq '    +os.path.join(options.childDir, 'inter', 'inter.outseq.rev') +\
                                     ' -chrname '  +chrom+     \
                                     ' -branchlength '  +str(options.stepSize)+ \
                                     ' -seed '     +options.seed+\
                                     ' -mes '      +os.path.join(options.childDir,'mobiles', 'mes.fa') + \
                                     ' -inannots ' +os.path.join(options.childDir, 'inter','inter.outannots.gff')+ \
                                     ' -statsfile '+os.path.join(options.childDir, 'stats',chrom+'.stats') +  \
                                     ' -codonsubs '+os.path.join(options.childDir, 'intra',  chrom+'.codonsubs') +\
                                     ' -outannots '+os.path.join(options.childDir, 'intra',  chrom+'.outannots.gff') + \
                                     ' -outgenome '+options.theChild + \
                                     ' -model '    +os.path.join(options.gParamsDir,'model.txt') + \
                                     ' -aln '      +os.path.join(options.childDir,  'intra',  chrom+'.aln.rev') + \
                                     ' -outseq '   +'LOCAL_DIR/'+chrom+'.outseq.rev' + \
                                     ' -log '      +os.path.join(options.childDir,  'logs', 'evo.'+chrom+'.log'))
        intraCMD +=LSC.commandPacker(CVT_BIN +\
                                     ' -fromrev ' +'LOCAL_DIR/'+chrom+'.outseq.rev' + \
                                     ' -tofasta ' +'LOCAL_DIR/'+chrom+'.outseq.fa' + \
                                     ' -log '     +os.path.join(options.childDir,  'logs', 'intra.'+chrom+'.tofasta.log'))
        intraCMD +=LSC.commandPacker(TRF_WRAP_BIN +\
                                     ' LOCAL_DIR/'+chrom+'.outseq.fa' + \
                                     ' 2 7 7 80 10 50 500 -d -h ')
        intraCMD +=LSC.commandPacker(ECHO_BIN +\
                                     ' \' \' ')
        intraCMD +=LSC.commandPacker(MV_BIN+\
                                     ' '+'LOCAL_DIR/'+chrom+'.*.dat '+ os.path.join(options.childDir, 'chr/'))
        intraCMD +=LSC.commandPacker(MV_BIN+\
                                     ' '+'LOCAL_DIR/'+chrom+'.outseq.rev '+ os.path.join(options.childDir, 'chr/')+' ')
        intraCMD +='"'
        newChild.attrib['command'] = intraCMD
    FILE.close()
    
    jobElm=xmlTree.getroot()
    followUpCommand = CYCLE_MAIN3 + \
                      ' --parent ' + options.parentDir +\
                      ' --child '  + options.childDir +\
                      ' --params ' + options.gParamsDir +\
                      ' --seed '   + options.seed+\
                      ' --jobFile JOB_FILE '
    jobElm.attrib['command'] = followUpCommand
    ##############################
    # finished
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)
