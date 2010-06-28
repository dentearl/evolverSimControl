#!/usr/bin/env python
"""
cycleMain_3.py
 dent earl, dearl (a) soe ucsc edu

cycleMain_3.py is a python wrapper for the evolver suite of genome
evolution tools. cycleMain_3 is the third in a series of four 
wrappers that is written to interface with jobTree.py, a cluster
interface written by Benedict Paten. Other members of the wrappers
are:
    cycleMain_1.py cycleMain_2.py cycleMain_3.py cycleMain_4.py

cycleMain_3.py Handles:
 1) Merge Steps.
*2) python trf2gff

* New features Feb 2010
19 October 2009
"""
########################################
import xml.etree.ElementTree as ET
import getopt, sys, os, subprocess
from optparse import OptionParser
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimCycle   as LSY

programs = ['evolver_evo', 'evolver_cvt', 'evolver_trf2gff.py',
            'simCtrl_cycleMain_4.py', 'cat', 'simCtrl_commandEval.py']
LSC.verifyPrograms(programs)
(EVO_BIN, CVT_BIN, TRF2GFF_BIN,
 CYCLE_END, CAT_BIN, CMD_EVAL_BIN) = programs

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
                         'main', 'cycleStep_2_cycleMain_2_end')
    LSC.subTypeTimestamp(os.path.join(options.childDir,'cycleInfo.xml'),
                         'main', 'cycleStep_3_cycleMain_3_start')
    
    xmlTree = ET.parse(options.jobFile)
    jobElm=xmlTree.getroot()
    childrenElm = xmlTree.find('children')

    FILE = open(os.path.join(options.childDir, 'inter', 'inter.chrnames.txt'))
    TRF_CMD      = CMD_EVAL_BIN+\
                   ' JOB_FILE "'+\
                   LSC.commandPacker(TRF2GFF_BIN +\
                   ' ' + os.path.join(options.childDir,'intra','*.dat')+\
                   ' > ' + os.path.join(options.childDir,'intra','trfannots.gff'))
    CAT_GFF_CMD  = CAT_BIN + ' ' + os.path.join(options.childDir,'intra','trfannots.gff')
    EVO_CMD  = ' '
    CVT_CMD  = ' '

    # build up the commands to make use of all the chromosome files that now exist.
    firstPass=True
    for line in FILE:
        chrom = line.rstrip()
        CAT_GFF_CMD = CAT_GFF_CMD +\
                      ' '+os.path.join(options.childDir, 'intra', chrom+'.outannots.gff')
        if firstPass: # these commands require comma separated values
            EVO_CMD = EVO_CMD+\
                     ' '+os.path.join(options.childDir, 'intra', chrom+'.aln.rev')
            CVT_CMD = CVT_CMD+\
                     ' '+os.path.join(options.childDir, 'intra', chrom+'.outseq.rev')
            firstPass=False
        else:
            EVO_CMD=EVO_CMD+\
                     ','+os.path.join(options.childDir, 'intra', chrom+'.aln.rev')
            CVT_CMD=CVT_CMD+\
                     ','+os.path.join(options.childDir, 'intra', chrom+'.outseq.rev')
    FILE.close()
    CAT_GFF_CMD= TRF_CMD +\
                 LSC.commandPacker(CAT_GFF_CMD + ' > '+os.path.join(options.childDir, 'intra', 'evannots.gff'))+'"'
    EVO_MERGE= CMD_EVAL_BIN +' JOB_FILE "'+\
               LSC.commandPacker(EVO_BIN+\
                                 ' -mergechrs '+EVO_CMD+ \
                                 ' -outgenome '+options.theChild+\
                                 ' -out '+os.path.join(options.childDir, 'intra', 'intra.aln.rev'))+'"'
    CVT_MERGE= CMD_EVAL_BIN +' JOB_FILE "'+\
               LSC.commandPacker(CVT_BIN+\
                                 ' -mergerevseqs '+CVT_CMD+\
                                 ' -out '+os.path.join(options.childDir, 'seq.rev'))+'"'

    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = CAT_GFF_CMD
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = EVO_MERGE
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = CVT_MERGE
    
    # follow up job will be the the final transalign step (cycleEnd.py)
    jobElm=xmlTree.getroot()
    followUpCommand = CYCLE_END +\
                      ' --parent ' + options.parentDir +\
                      ' --child '  + options.childDir +\
                      ' --params ' + options.gParamsDir +\
                      ' --seed '   + options.seed+\
                      ' --jobFile JOB_FILE'
    jobElm.attrib['command'] = followUpCommand
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)
