#!/usr/bin/env python
"""
simCtrl_stats_2.py
dent earl, dearl (a) soe ucsc edu
11 mar 2010

This script performs a series of statistical tests
on an evolver cycle. The script is written such that
it can be run either as a part of an evolver cycle,
tacked onto the end, or at a later time in some
separate process.

"""
##############################
import xml.etree.ElementTree as ET
from optparse import OptionParser
import os, sys
#import eval.lib.libSimTree as LST
import eval.lib.libSimControl as LSC
import eval.lib.libSimStats as LSS

programs = ['simCtrl_commandEval.py', 'evolver_evo', 'evolver_gff_cdsutr2exons.py',
            'evolver_gff_exons2introns.py', 'evolver_gff_featurestats2.sh',
            'cat', 'egrep', 'simCtrl_cycleStats_4.py', 'evolver_cvt',
            'evolver_gff_featurestats2.py','evolver_codon_report.pl', 'evolver_drawrev',
            'evolver_merge_evostats.py', 'evolver_evostats_report.py', 'touch',
            'simCtrl_completeTimestamp.py']
LSC.verifyPrograms(programs)
(CMD_EVAL_BIN, EVO_BIN, CDSUTR2EXON_BIN, EXON2INTRON_BIN,
 FEATSTATS_SH, CAT_BIN, EGREP_BIN, STAT_4_BIN, CVT_BIN,
 FEATSTATS_PY, CODONREPORT_BIN, DRAWREV_BIN, MERGE_BIN,
 EVOSTATS_BIN, TOUCH, TIMESTAMP_BIN)  = programs

DRAW_REV_BLOCK_SIZE=10000
DRAW_REV_NT_PER_PIX=10000

def usage():
    print 'USAGE: %s --childDir <dir> --parentDir <dir> --jobFile <JOB_FILE> ' %(sys.argv[0])
    print __doc__
    sys.exit(2)

def main():
    parser=OptionParser()
    LSS.initOptions(parser)
    (options, args) = parser.parse_args()
    LSS.checkOptions(options)
    LSC.subTypeTimestamp(os.path.join(options.childDir, 'cycleInfo.xml'), 'stats', 'cycleStep_7_cycleStats_3_end')
    LSC.subTypeTimestamp(os.path.join(options.childDir, 'cycleInfo.xml'), 'stats', 'cycleStep_8_cycleStats_4_start')
    
    xmlTree = ET.parse(options.jobFile)
    jobElm=xmlTree.getroot()
    childrenElm = xmlTree.find('children')

    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(CAT_BIN+\
                                ' '+ os.path.join(options.childDir, 'stats', 'tmpstats.root.difflength.tmp')+\
                                ' '+ os.path.join(options.childDir, 'stats', 'tmpstats.root.diffcompost.tmp')+\
                                ' '+ os.path.join(options.childDir, 'stats', 'tmpstats.root.diffannots.tmp')+\
                                ' > '+ os.path.join(options.childDir, 'stats', 'childRoot.diff.txt'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD
    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(MERGE_BIN+\
                                ' '+os.path.join(options.childDir,'stats','merged_cycle.stats')+\
                                ' '+os.path.join(options.parentDir,'stats','merged_root.stats')+\
                                ' > '+os.path.join(options.childDir,'stats','merged_root.stats'))+\
              LSC.commandPacker(EVOSTATS_BIN+\
                                ' '+os.path.join(options.childDir,'stats','merged_root.stats')+\
                                ' > '+os.path.join(options.childDir,'stats','root_events.txt'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD

    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(EVO_BIN+\
                                ' -cdsalns '+ os.path.join(options.childDir, 'root.aln.rev')+ \
                                ' -alns '   + os.path.join(options.childDir, 'stats', 'cds_alns.rev')+ \
                                ' -annots1 '+ os.path.join(options.rootDir, 'stats', 'cds_annots.gff')+ \
                                ' -annots2 '+ os.path.join(options.childDir, 'annots.gff')+ \
                                ' -log '+ os.path.join(options.childDir, 'logs', 'root.cds_alns.log'))+\
              LSC.commandPacker(EVO_BIN+\
                                ' -getcodonsubs '+os.path.join(options.childDir,'stats', 'cds_alns.rev')+\
                                ' -out '+os.path.join(options.childDir, 'stats', 'root.codon_sub.txt')+\
                                ' -log '+os.path.join(options.childDir, 'logs','root.getcodonsubs.log'))+\
              LSC.commandPacker(CODONREPORT_BIN+\
                                ' '+os.path.basename(options.rootDir)+\
                                ' '+os.path.basename(options.childDir)+\
                                ' < '+os.path.join(options.childDir, 'stats', 'root.codon_sub.txt')+\
                                ' > '+os.path.join(options.childDir, 'stats', 'root.protstats.txt'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD
    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(DRAWREV_BIN+\
                                ' -fromrev '+os.path.join(options.childDir, 'root.aln.rev')+\
                                ' -tocmap '+os.path.join(options.childDir, 'stats', 'img.root.cmap.pdf')+\
                                ' -blocksize '+str(DRAW_REV_BLOCK_SIZE))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD
    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(DRAWREV_BIN+\
                                ' -fromrev '+os.path.join(options.childDir, 'root.aln.rev')+\
                                ' -tolmap ' +os.path.join(options.childDir, 'stats', 'img.root.lmap.png')+\
                                ' -npp '    +str(DRAW_REV_NT_PER_PIX)) +'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD

    followUpCommand = CMD_EVAL_BIN+\
                      ' JOB_FILE "'+\
                      LSC.commandPacker(TIMESTAMP_BIN+\
                                        ' --cycleDir '+options.childDir+\
                                        ' --timeType=stats ')+'"'
    jobElm.attrib['command'] = followUpCommand
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main()
