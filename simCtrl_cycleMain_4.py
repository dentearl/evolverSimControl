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
########################################
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

# def transCMDBuilder(options):
#     transCMD = CMD_EVAL_BIN+\
#                ' --statXML '+os.path.join(options.childDir, 'logs','trans.info.xml')+\
#                ' JOB_FILE "'+\
#                LSC.commandPacker(TRANS_BIN +\
#                                  ' -in1 '+os.path.join(options.childDir, 'inter','inter.aln.rev')+ \
#                                  ' -in2 '+os.path.join(options.childDir, 'intra', 'intra.aln.rev')+ \
#                                  ' -out '+os.path.join(options.childDir, 'inter-intra.aln.rev')+ \
#                                  ' -log '+os.path.join(options.childDir, 'logs', 'transalign1.log'))
#     if( os.path.isfile(os.path.join(options.parentDir, 'root.aln.rev'))):
#         transCMD += LSC.commandPacker(TRANS_BIN +\
#                                      ' -in1 '+os.path.join(options.parentDir,'root.aln.rev')+ \
#                                      ' -in2 '+os.path.join(options.childDir, 'inter-intra.aln.rev')+ \
#                                      ' -out '+os.path.join(options.childDir, 'root.aln.rev')+ \
#                                      ' -log '+os.path.join(options.childDir, 'logs', 'transalign2.log'))
#     else:
#         # base case, the parent *is* the root.
#         transCMD += LSC.commandPacker(LINK_BIN +\
#                                       ' -s '+os.path.join(options.childDir,'inter-intra.aln.rev')+\
#                                       ' '+os.path.join(options.childDir,'root.aln.rev'))
#     transCMD +='"'
#     return transCMD

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
    # Transalign ONLY if this is a LEAF genome,
    # i.e., it is the final genome in a simulation.
    ########################################
#     if options.isLeaf:
#         transCMD = transCMDBuilder(options)
#         newChild = ET.SubElement(childrenElm, 'child')
#         newChild.attrib['command'] = transCMD

    ####################
    # gene decativation step
    gDActCMD = CMD_EVAL_BIN+\
               ' JOB_FILE "'+\
               LSC.commandPacker(GDACT_BIN +\
                                 ' ' + os.path.join(options.parentDir, 'annots.gff') +\
                                 ' ' + os.path.join(options.childDir, 'intra', 'evannots.gff')+\
                                 ' ' + os.path.join(options.childDir, 'annots.gff')+\
                                 ' ' + EVO_BIN+\
                                 ' >& ' + os.path.join(options.childDir, 'logs', 'gene_deactivate.log'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = gDActCMD

    # Last step, add the cycle timestamp.
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
