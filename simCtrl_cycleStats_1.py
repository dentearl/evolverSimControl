#!/usr/bin/env python
"""
simCtrl_stats_1.py
dent earl, dearl (a) soe ucsc edu
9 mar 2010

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
            'cat', 'egrep', 'simCtrl_cycleStats_2.py', 'evolver_cvt',
            'evolver_gff_featurestats2.py']
LSC.verifyPrograms(programs)
(CMD_EVAL_BIN, EVO_BIN, CDSUTR2EXON_BIN, EXON2INTRON_BIN,
 FEATSTATS_SH, CAT_BIN, EGREP_BIN, STAT_2_BIN, CVT_BIN, FEATSTATS_PY)  = programs

def usage():
    print 'USAGE: %s --cycleDir <dir> --parentDir <dir> --jobFile <JOB_FILE> ' %(sys.argv[0])
    print __doc__
    sys.exit(2)

def main(argv):
    parser=OptionParser()
    LSS.initOptions( parser )
    ( options, args ) = parser.parse_args()
    LSS.checkOptions( options, usage )
    LSC.typeTimestamp(os.path.join(options.childDir,'cycleInfo.xml'),
                      'stats', 'Start')
    LSC.subTypeTimestamp(os.path.join(options.childDir,'cycleInfo.xml'),
                         'stats', 'cycleStep_5_cycleStats_1_start')
    
    xmlTree = ET.parse(options.jobFile)
    jobElm=xmlTree.getroot()
    childrenElm = xmlTree.find('children')
    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(EVO_BIN+\
                                ' -probstats '+os.path.join(options.childDir, 'annots.gff')+\
                                ' -log '+os.path.join(options.childDir,'stats','probstats.txt'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD
    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(EVO_BIN+\
                                ' -annotstats '+os.path.join(options.childDir, 'annots.gff')+\
                                ' -seq '+os.path.join(options.childDir, 'seq.rev')+\
                                ' -log '+os.path.join(options.childDir, 'stats','annotstats.txt'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD
    statCMD = CMD_EVAL_BIN+\
              ' JOB_FILE "'+\
              LSC.commandPacker(EGREP_BIN+\
                                ' \'CDS|UTR\' '+os.path.join(options.childDir,'annots.gff')+\
                                ' > '+os.path.join(options.childDir,'stats','cds_annots.gff'))+\
              LSC.commandPacker(CDSUTR2EXON_BIN+\
                                ' '+os.path.join(options.childDir,'stats','cds_annots.gff')+\
                                ' > '+os.path.join(options.childDir,'stats','exons.gff'))+\
              LSC.commandPacker(EXON2INTRON_BIN+\
                                ' '+os.path.join(options.childDir,'stats','exons.gff')+\
                                ' > '+os.path.join(options.childDir,'stats','introns.gff'))+\
              LSC.commandPacker(CAT_BIN+\
                                ' '+os.path.join(options.childDir,'annots.gff')+\
                                ' '+os.path.join(options.childDir,'stats','exons.gff')+\
                                ' '+os.path.join(options.childDir,'stats','introns.gff')+\
                                ' > '+os.path.join(options.childDir,'stats','expanded_annots.gff'))+\
              LSC.commandPacker(FEATSTATS_SH+\
                                ' '+FEATSTATS_PY+\
                                ' '+os.path.join(options.parentDir,'stats','expanded_annots.gff')+\
                                ' '+os.path.join(options.childDir,'stats','expanded_annots.gff')+\
                                ' '+os.path.basename(options.parentDir)+\
                                ' '+os.path.basename(options.childDir)+\
                                ' '+os.path.join(options.parentDir,'seq.rev')+\
                                ' '+os.path.join(options.childDir,'seq.rev')+\
                                ' '+CVT_BIN+\
                                ' > '+os.path.join(options.childDir,'stats','tmpstats.diffannots.tmp'))+'"'
    newChild = ET.SubElement(childrenElm, 'child')
    newChild.attrib['command'] = statCMD

    followUpCommand  = STAT_2_BIN
    followUpCommand += ' --childDir ' + options.childDir
    followUpCommand += ' --parentDir '+ options.parentDir
    if options.noMEs:
        followUpCommand += ' --noMEs'
    followUpCommand += ' --jobFile JOB_FILE '
    jobElm.attrib['command'] = followUpCommand
    xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)
