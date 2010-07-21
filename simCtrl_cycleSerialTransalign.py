#!/usr/bin/env python
"""
cycleSerialTransalign.py
 dent earl, dearl (a) soe ucsc edu

cycleSerialTransalign.py is a python wrapper for the evolver suite of
genome evolution tools. This script is executed for a given parent
and child cycle pair, but it is executed in parallel with the grandchild
cycle. Thus, it is transaligning the previous cycle (the parent) with it's
parent (the grandparent). This is done to increase parallelization.

cycleSerialTransalign.py Handles:
 1) 1st Transalign step.
 2) 2nd Transalign step.

19 July 2010
"""
########################################
import xml.etree.ElementTree as ET
import os, subprocess, sys
from optparse import OptionParser
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimCycle   as LSY

programs = ['evolver_transalign', 'touch', 'simCtrl_commandEval.py', 'evolver_gene_deactivate.sh',
            'simCtrl_cycleStats_1.py', 'simCtrl_completeTimestamp.py', 'ln', 'evolver_evo',
            'evolver_codon_report.pl', 'cat', 'evolver_drawrev']
LSC.verifyPrograms(programs)
(TRANS_BIN, TOUCH, CMD_EVAL_BIN, GDACT_BIN,
 STATS_BIN, TIMESTAMP_BIN, LINK_BIN, EVO_BIN,
 CODONREPORT_BIN, CAT_BIN, DRAWREV_BIN) = programs

DRAW_REV_BLOCK_SIZE=10000
DRAW_REV_NT_PER_PIX=100000

def usage():
    print "USAGE: %s --parent parentDir/ --child childDir --jobFile JOB_FILE [optional: --step ]" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def findGrandParent(options):
    # open the parent dir's cycleInfo file
    file = os.path.join(options.parentDir, 'cycleInfo.xml')
    infoTree = ET.parse(file)
    root = infoTree.getroot()
    # find the parent of the parent, record it as the grandParent
    tObj = root.find('parentDir')
    options.grandParentDir = str(tObj.text)

def customOptions(parser):
    parser.add_option('--isLeaf', action='store_true', dest='isLeaf',
                      default=False, help='Is this a leaf run? If so, do not talign grand parent -> parent, talign parent -> child')

def main(argv):
    parser=OptionParser()
    LSC.standardOptions(parser)
    LSY.standardOptions(parser)
    customOptions(parser)
    (options, args) = parser.parse_args()
    LSC.standardOptionsCheck(options, usage)
    LSY.standardOptionsCheck(options, usage)
    findGrandParent(options)
#    LSC.subTypeTimestamp(os.path.join(options.parentDir,'cycleInfo.xml'),
#                         'talign', 'cycleStep_3_cycleMain_3_end')
#    LSC.subTypeTimestamp(os.path.join(options.parentDir,'cycleInfo.xml'),
#                         'talign', 'cycleStep_4_cycleMain_4_start')

    xmlTree = ET.parse(options.jobFile)
    jobElm=xmlTree.getroot()
    childrenElm = xmlTree.find('children')

    if not options.isLeaf:
        A = options.grandParentDir
        B = options.parentDir
    else:
        A = options.parentDir
        B = options.childDir
    
    ########################################
    # Transalign
    ########################################
    # the two transalign commands will be performed in series.
    transCMD = CMD_EVAL_BIN+\
               ' --statXML '+os.path.join(B, 'logs','trans.info.xml')+\
               ' JOB_FILE "'
    # command_0
    transCMD +=LSC.commandPacker(TRANS_BIN +\
                                 ' -in1 '+os.path.join(B, 'inter','inter.aln.rev')+ \
                                 ' -in2 '+os.path.join(B, 'intra', 'intra.aln.rev')+ \
                                 ' -out '+os.path.join(B, 'inter-intra.aln.rev')+ \
                                 ' -log '+os.path.join(B, 'logs', 'transalign1.log'))
    # command_1
    if( os.path.isfile(os.path.join(A, 'root.aln.rev'))):
        transCMD += LSC.commandPacker(TRANS_BIN +\
                                     ' -in1 '+os.path.join(A,'root.aln.rev')+ \
                                     ' -in2 '+os.path.join(B, 'inter-intra.aln.rev')+ \
                                     ' -out '+os.path.join(B, 'root.aln.rev')+ \
                                     ' -log '+os.path.join(B, 'logs', 'transalign2.log'))
    else:
        # base case, the parent *is* the root.
        transCMD += LSC.commandPacker(LINK_BIN +\
                                      ' -s '+os.path.join(B,'inter-intra.aln.rev')+\
                                      ' '+os.path.join(B,'root.aln.rev'))
    # command_2
    transCMD += LSC.commandPacker(EVO_BIN+\
                                  ' -cdsalns '+ os.path.join(B, 'intra', 'intra.aln.rev')+ \
                                  ' -alns '   + os.path.join(B, 'stats', 'cds_alns.rev')+ \
                                  ' -annots1 '+ os.path.join(B, 'inter', 'inter.outannots.gff')+ \
                                  ' -annots2 '+ os.path.join(B, 'annots.gff')+ \
                                  ' -log '+ os.path.join(B, 'logs', 'cds_alns.log'))
    # command_3
    transCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                 ' -fromrev '+os.path.join(B, 'inter', 'inter.aln.rev')+\
                                 ' -tocmap '+os.path.join(B, 'stats', 'img.cycle.cmap.pdf')+\
                                 ' -blocksize '+str(DRAW_REV_BLOCK_SIZE))
    # command_4
    transCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                 ' -fromrev '+os.path.join(B, 'inter', 'inter.aln.rev')+\
                                 ' -tolmap ' +os.path.join(B, 'stats', 'img.cycle.lmap.png')+\
                                 ' -npp '    +str(DRAW_REV_NT_PER_PIX))
    # command_5
    transCMD +=LSC.commandPacker(EVO_BIN+\
                                 ' -nologcmdlineandtime '+\
                                 ' -ancstats '+os.path.join(B, 'root.aln.rev') +\
                                 ' -log '+ os.path.join(B, 'stats', 'tmpstats.root.difflength.tmp'))
    # command_6
    transCMD +=LSC.commandPacker(EVO_BIN+\
                                 ' -nologcmdlineandtime '+\
                                 ' -ancstats '+os.path.join(B, 'intra', 'intra.aln.rev') +\
                                 ' -log '+ os.path.join(B, 'stats', 'tmpstats.difflength.tmp'))
    # command_7
    transCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                 ' -fromrev '+os.path.join(B, 'root.aln.rev') +\
                                 ' -tocmap '+os.path.join(B, 'stats', 'img.root.cmap.pdf') +\
                                 ' -blocksize '+str(DRAW_REV_BLOCK_SIZE))
    # command_8
    transCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                 ' -fromrev '+os.path.join(B, 'root.aln.rev')+\
                                 ' -tolmap ' +os.path.join(B, 'stats', 'img.root.lmap.png')+\
                                 ' -npp '    +str(DRAW_REV_NT_PER_PIX))
    # command_9
    transCMD +=LSC.commandPacker(EVO_BIN+\
                                 ' -cdsalns '+ os.path.join(B, 'root.aln.rev')+ \
                                 ' -alns '   + os.path.join(B, 'stats', 'cds_alns.rev')+ \
                                 ' -annots1 '+ os.path.join(options.rootDir, 'stats', 'cds_annots.gff')+ \
                                 ' -annots2 '+ os.path.join(B, 'annots.gff')+ \
                                 ' -log '+ os.path.join(B, 'logs', 'root.cds_alns.log'))
    # command_10
    transCMD +=LSC.commandPacker(EVO_BIN+\
                                 ' -getcodonsubs '+os.path.join(B,'stats', 'cds_alns.rev')+\
                                 ' -out '+os.path.join(B, 'stats', 'codon_sub.txt')+\
                                 ' -log '+os.path.join(B, 'logs','getcodonsubs.log'))
    # command_11
    transCMD +=LSC.commandPacker(EVO_BIN+\
                                 ' -getcodonsubs '+os.path.join(B,'stats', 'cds_alns.rev')+\
                                 ' -out '+os.path.join(B, 'stats', 'root.codon_sub.txt')+\
                                 ' -log '+os.path.join(B, 'logs','root.getcodonsubs.log'))
    # command_12
    transCMD +=LSC.commandPacker(CODONREPORT_BIN+\
                                 ' '+os.path.basename(A)+\
                                 ' '+os.path.basename(B)+\
                                 ' < '+os.path.join(B, 'stats', 'codon_sub.txt')+\
                                 ' > '+os.path.join(B, 'stats', 'protstats.txt'))
    # command_13
    transCMD +=LSC.commandPacker(CODONREPORT_BIN+\
                                 ' '+os.path.basename(options.rootDir)+\
                                 ' '+os.path.basename(B)+\
                                 ' < '+os.path.join(B, 'stats', 'root.codon_sub.txt')+\
                                 ' > '+os.path.join(B, 'stats', 'root.protstats.txt'))
    # command_14
    transCMD +=LSC.commandPacker(CAT_BIN+\
                                 ' '+ os.path.join(B, 'stats', 'tmpstats.difflength.tmp')+\
                                 ' '+ os.path.join(B, 'stats', 'tmpstats.diffcompost.tmp')+\
                                 ' '+ os.path.join(B,'stats','tmpstats.diffannots.tmp')+\
                                 ' > '+ os.path.join(B, 'stats', 'childParent.diff.txt'))
    # command_15
    transCMD +=LSC.commandPacker(CAT_BIN+\
                                ' '+ os.path.join(B, 'stats', 'tmpstats.root.difflength.tmp')+\
                                ' '+ os.path.join(B, 'stats', 'tmpstats.root.diffcompost.tmp')+\
                                ' '+ os.path.join(B, 'stats', 'tmpstats.root.diffannots.tmp')+\
                                ' > '+ os.path.join(B, 'stats', 'childRoot.diff.txt'))
    transCMD +='"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = transCMD

    # Last step, add the cycle timestamp.
#     followUpCommand = CMD_EVAL_BIN+\
#                       ' JOB_FILE "'+\
#                       LSC.commandPacker(TIMESTAMP_BIN+\
#                                         ' --cycleDir '+options.parentDir+\
#                                         ' --timeType=talign '+\
#                       LSC.commandPacker(STATS_BIN+\
#                                         ' --childDir '+options.parentDir+\
#                                         ' --parentDir '+options.grandParentDir+\
#                                         ' --jobFile JOB_XML'))+'"'
#    jobElm=xmlTree.getroot()
#    jobElm.attrib['command'] = followUpCommand
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)


