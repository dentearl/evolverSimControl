#!/usr/bin/env python
"""
cycleMain_4.py
 dent earl, dearl (a) soe ucsc edu

cycleMain_4.py is a python wrapper for the evolver suite of genome
evolution tools. cycleMain_4 is the fourth in a series of four
wrappers that is written to interface with jobTree.py, a cluster
interface written by Benedict Paten. Other members of the wrappers
are:
    cycleMain_1.py cycleMain_2.py cycleMain_3.py cycleMain_4.py

cycleMain_4.py Handles:
 1) 2nd Transalign step.
 2) CDS align step.

7 October 2009
"""
########################################
import xml.etree.ElementTree as ET
import os, subprocess, sys
from optparse import OptionParser
import eval.lib.libSimControl as LSC
import eval.lib.libSimCycle   as LSY

programs = ['evolver_evo', 'evolver_cvt', 'evolver_transalign',
            'touch', 'simCtrl_commandEval.py', 'evolver_gene_deactivate.sh',
            'simCtrl_cycleStats_1.py', 'simCtrl_completeTimestamp.py']
LSC.verifyPrograms(programs)
(EVO_BIN, CVT_BIN, TRANS_BIN, TOUCH, CMD_EVAL_BIN, GDACT_BIN,
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

    ########################################
    # Transalign 2, CDSalign PARALLEL 3
    ########################################
    transCMD = CMD_EVAL_BIN+\
               ' JOB_FILE "'+\
               LSC.commandPacker(TRANS_BIN +\
                                 ' -in1 '+os.path.join(options.childDir, 'inter', options.theParent+'.inter.aln.rev')+ \
                                 ' -in2 '+os.path.join(options.childDir, 'chr', 'intra.aln.rev')+ \
                                 ' -out '+os.path.join(options.childDir, 'root.aln.rev')+ \
                                 ' -log '+os.path.join(options.childDir, 'logs', 'transalign.log'))+'"'
    gDActCMD = CMD_EVAL_BIN+\
               ' JOB_FILE "'+\
               LSC.commandPacker(GDACT_BIN +\
                                 ' ' + os.path.join(options.parentDir, 'annots.gff') +\
                                 ' ' + os.path.join(options.childDir, 'chr', 'evannots.gff')+\
                                 ' ' + os.path.join(options.childDir, 'annots.gff')+\
                                 ' ' + EVO_BIN+\
                                 ' >& ' + os.path.join(options.childDir, 'logs', 'gene_deactivate.log'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = transCMD
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = gDActCMD

    # Last step, add the timestamp.
    followUpCommand = CMD_EVAL_BIN+\
                      ' JOB_FILE "'+\
                      LSC.commandPacker(TIMESTAMP_BIN+\
                                        ' --cycleDir '+options.childDir+\
                                        ' --timeType=main '+\
                      LSC.commandPacker(STATS_BIN+\
                                        ' --childDir '+options.childDir+\
                                        ' --parentDir '+options.parentDir+\
                                        ' --jobFile JOB_XML'))+'"'
                                        
    jobElm=xmlTree.getroot()
    jobElm.attrib['command'] = followUpCommand
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)
