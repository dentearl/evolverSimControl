#!/usr/bin/env python
"""
simCtrl_stats_4.py
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
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimStats as LSS

programs = ['simCtrl_commandEval.py', 'evolver_evo', 'evolver_gff_cdsutr2exons.py',
            'evolver_gff_exons2introns.py', 'evolver_gff_featurestats2.sh',
            'cat', 'egrep', 'simCtrl_cycleStats_4.py', 'evolver_cvt',
            'evolver_gff_featurestats2.py',
            'evolver_merge_evostats.py', 'evolver_evostats_report.py', 'touch',
            'simCtrl_completeTimestamp.py']
LSC.verifyPrograms(programs)
(CMD_EVAL_BIN, EVO_BIN, CDSUTR2EXON_BIN, EXON2INTRON_BIN,
 FEATSTATS_SH, CAT_BIN, EGREP_BIN, STAT_4_BIN, CVT_BIN,
 FEATSTATS_PY, MERGE_BIN,
 EVOSTATS_BIN, TOUCH, TIMESTAMP_BIN)  = programs

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
              ' JOB_FILE "'
    statCMD +=LSC.commandPacker(MERGE_BIN+\
                                ' '+os.path.join(options.childDir,'stats','merged_cycle.stats')+\
                                ' '+os.path.join(options.parentDir,'stats','merged_root.stats')+\
                                ' > '+os.path.join(options.childDir,'stats','merged_root.stats'))
    statCMD +=LSC.commandPacker(EVOSTATS_BIN+\
                                ' '+os.path.join(options.childDir,'stats','merged_root.stats')+\
                                ' > '+os.path.join(options.childDir,'stats','root_events.txt'))
    statCMD +='"'
    
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
