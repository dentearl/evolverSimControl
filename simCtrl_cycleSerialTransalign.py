#!/usr/bin/env python
"""
cycleSerialTransalign.py
 dent earl, dearl (a) soe ucsc edu

cycleSerialTransalign.py is a python wrapper for the evolver suite of
genome evolution tools. This script is executed for a given target
cycle.

cycleSerialTransalign.py Handles:
 1) 1st Transalign step.
 2) 2nd Transalign step.
 3) cdsaln intra to inter
 4) drawrev step 1
 5) drawrev step 2
 6) anc stats step root
 7) anc stats step intra
 8) cdsaln to root
 9) codon sub stats cdsaln inter-intra
10) codon sub stats cdsaln root
11) codon report inter-intra parsing
12) codon report root parsing
13) cat together within cycle stats
14) cat together between cycle stats

19 July 2010
"""
########################################
import xml.etree.ElementTree as ET
from sonLib.bioio import newickTreeParser
import os, subprocess, sys
from optparse import OptionParser
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

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
    print "USAGE: %s --targetDir cycle/ --jobFile JOB_FILE" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def findParent(options):
    # open the parent dir's cycleInfo file
    file = os.path.join(options.targetDir, 'cycleInfo.xml')
    infoTree = ET.parse(file)
    root = infoTree.getroot()
    # find the parent of the parent, record it as the grandParent
    tObj = root.find('parentDir')
    options.parentDir = os.path.abspath( str(tObj.text) )

def parentWasBranchOrRoot( options ):
    """ returns True if the parent cycle was either a branch point or the root
    """
    if os.path.samefile( options.parentDir, options.rootDir ):
        return True
    file = os.path.join(options.parentDir, 'cycleInfo.xml')
    infoTree = ET.parse(file)
    root = infoTree.getroot()
    if len( root.findall('followUpCommand') ) == 3:
        # 3 is a branch, 2 is a stem
        return True
    return False

def customOptions(parser):
    parser.add_option('-t', '--targetDir',dest='targetDir',
                      help='Cycle directory to perform the transalign upon.')
    parser.add_option('-j', '--jobFile',dest='jobFile',
                      help='jobFile, passed in by jobTree.py.')
    parser.add_option('--debug', action='store_true', dest='isDebug',
                      default=False, help='Performs no operations, prints out status of variables.')
    
        
def customOptionsCheck( options, usage ):
    if options.targetDir == None:
        sys.stderr.write('%s: Error, specify --targetDir.\n' % sys.argv[0])
        usage()
    if not os.path.isdir(options.targetDir):
        sys.stderr.write('%s: Error, Parent dir "%s" not a directory!\n' % (sys.argv[0], options.targetDir))
        usage()
    options.targetDir  = os.path.abspath(options.targetDir)
    (options.simDir, tail) = os.path.split(options.targetDir)

    xmlTree = ET.parse( os.path.join( options.simDir, 'simulationInfo.xml' ) )
    tree = xmlTree.find('tree')
    nt = newickTreeParser( tree.text, 0.0 )
    # rootName = LST.newickRootName( nt )
    rootNameObj = xmlTree.find('rootDir')
    options.rootName = os.path.basename( rootNameObj.text )
    options.rootDir    = os.path.abspath( rootNameObj.text )
    if options.jobFile == None:
        sys.stderr.write('%s: Error, specify --jobFile.\n' % sys.argv[0])
        usage()

def main(argv):
    parser=OptionParser()
    customOptions(parser)
    (options, args) = parser.parse_args()
    if not options.isDebug:
        customOptionsCheck(options, usage)
    else:
        options.rootDir = 'root'
    findParent(options)

    A = options.parentDir
    B = options.targetDir
    
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
    if parentWasBranchOrRoot( options ):
        # In these cases the alignment above the branch point should not be carried
        # into the the descendant genomes. Alignments should only go back to the most
        # recent branch point.
        transCMD += LSC.commandPacker(LINK_BIN +\
                                      ' -s '+os.path.join(B,'inter-intra.aln.rev')+\
                                      ' '+os.path.join(B,'aln.rev'))
    else:
        transCMD += LSC.commandPacker(TRANS_BIN +\
                                      ' -in1 '+os.path.join(A, 'aln.rev')+ \
                                      ' -in2 '+os.path.join(B, 'inter-intra.aln.rev')+ \
                                      ' -out '+os.path.join(B, 'aln.rev')+ \
                                      ' -log '+os.path.join(B, 'logs', 'transalign2.log'))
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
                                 ' -ancstats '+os.path.join(B, 'aln.rev') +\
                                 ' -log '+ os.path.join(B, 'stats', 'tmpstats.root.difflength.tmp'))
    # command_6
    transCMD +=LSC.commandPacker(EVO_BIN+\
                                 ' -nologcmdlineandtime '+\
                                 ' -ancstats '+os.path.join(B, 'intra', 'intra.aln.rev') +\
                                 ' -log '+ os.path.join(B, 'stats', 'tmpstats.difflength.tmp'))
    # command_7
    transCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                 ' -fromrev '+os.path.join(B, 'aln.rev') +\
                                 ' -tocmap '+os.path.join(B, 'stats', 'img.root.cmap.pdf') +\
                                 ' -blocksize '+str(DRAW_REV_BLOCK_SIZE))
    # command_8
    transCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                 ' -fromrev '+os.path.join(B, 'aln.rev')+\
                                 ' -tolmap ' +os.path.join(B, 'stats', 'img.root.lmap.png')+\
                                 ' -npp '    +str(DRAW_REV_NT_PER_PIX))
    # command_9
    transCMD +=LSC.commandPacker(EVO_BIN+\
                                 ' -cdsalns '+ os.path.join(B, 'aln.rev')+ \
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
    if options.isDebug:
        print '#'*20 + '\ntransCMD is: \n%s\n' % '\n'.join(LSC.commandUnPacker( transCMD ))
    else:
        xmlTree = ET.parse(options.jobFile)
        jobElm=xmlTree.getroot()
        childrenElm = xmlTree.find('children')
        newChild = ET.SubElement(childrenElm, 'child')
        newChild.attrib['command'] = transCMD
        xmlTree.write(options.jobFile)

if __name__ == "__main__":
    main(sys.argv)


