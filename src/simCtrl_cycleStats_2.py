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
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimStats as LSS

programs = ['simCtrl_commandEval.py', 'evolver_evo', 'evolver_gff_cdsutr2exons.py',
            'evolver_gff_exons2introns.py', 'evolver_gff_featurestats2.sh',
            'cat', 'egrep', 'simCtrl_cycleStats_3.py', 'evolver_cvt',
            'evolver_gff_featurestats2.py',
            'evolver_merge_evostats.py', 'evolver_evostats_report.py',
            'evolver_mobile_report.pl']
LSC.verifyPrograms(programs)
(CMD_EVAL_BIN, EVO_BIN, CDSUTR2EXON_BIN, EXON2INTRON_BIN,
 FEATSTATS_SH, CAT_BIN, EGREP_BIN, STAT_3_BIN, CVT_BIN,
 FEATSTATS_PY, MERGE_BIN,
 EVOSTATS_BIN, MOBILE_REPORT_BIN)  = programs

def usage():
    print 'USAGE: %s --childDir <dir> --parentDir <dir> --jobFile <JOB_FILE> ' %(sys.argv[0])
    print __doc__
    sys.exit(2)

def main():
    parser=OptionParser()
    LSS.initOptions(parser)
    (options, args) = parser.parse_args()
    LSS.checkOptions( options, usage )
    LSC.subTypeTimestamp(os.path.join(options.childDir, 'cycleInfo.xml'), 'stats', 'cycleStep_5_cycleStats_1_end')
    LSC.subTypeTimestamp(os.path.join(options.childDir, 'cycleInfo.xml'), 'stats', 'cycleStep_6_cycleStats_2_start')
    
    xmlTree = ET.parse(options.jobFile)
    jobElm=xmlTree.getroot()
    childrenElm = xmlTree.find('children')

    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(MERGE_BIN+\
                                ' '+os.path.join(options.childDir,'stats','*stats')+\
                                ' > '+os.path.join(options.childDir,'stats','merged_cycle.stats'))+\
              LSC.commandPacker(EVOSTATS_BIN+\
                                ' '+os.path.join(options.childDir,'stats','merged_cycle.stats')+\
                                ' > '+os.path.join(options.childDir,'stats','events.txt'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD
    if not options.noMEs:
        statCMD = CMD_EVAL_BIN+\
                  ' JOB_FILE "'+\
                  LSC.commandPacker(MOBILE_REPORT_BIN+\
                                    ' '+os.path.join(options.childDir,'logs','intra.*.log')+\
                                    ' --mobilesLog '+os.path.join(options.childDir,'logs','mobiles.log')+\
                                    ' > '+os.path.join(options.childDir,'stats','stats.mobiles.txt'))+'"'
        newChild = ET.SubElement(childrenElm, 'child')
        newChild.attrib['command'] = statCMD
        
    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(EVO_BIN+\
                                ' -nologcmdlineandtime '+\
                                ' -compost1 '+os.path.join(options.parentDir, 'seq.rev')+ \
                                ' -compost2 '+os.path.join(options.childDir, 'seq.rev')+ \
                                ' -log '+ os.path.join(options.childDir, 'stats', 'tmpstats.diffcompost.tmp'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD

    followUpCommand  = STAT_3_BIN
    followUpCommand += ' --childDir ' + options.childDir
    followUpCommand += ' --parentDir '+ options.parentDir
    if options.noMEs:
        followUpCommand += ' --noMEs'
    followUpCommand += ' --jobFile JOB_FILE '
    jobElm.attrib['command'] = followUpCommand
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main()