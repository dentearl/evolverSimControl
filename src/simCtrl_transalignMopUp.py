#!/usr/bin/env python
"""
transalignMopUp.py
 dent earl, dearl (a) soe ucsc edu

transalignMopUp.py is a python script that is meant to be run
in a jobTree command. Its purpose is to recover failed
transalign jobs (which were launched by cycleSerialTransalign.py).
This script is needed because there are instances during some
simulations where the cycleSerialTransalign.py job from
one cycle takes long enough that the cycleSerialTransalign for
the next cycle starts and fails. As soon as one of these fails it
cascades through the simulation.
The fix for this bug is going to be a complicated design fix,
so in the meantime I have this mopup script that you are reading.

13 December 2010
"""
##################################################
# Copyright (C) 2009-2011 by
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##################################################
import xml.etree.ElementTree as ET
from sonLib.bioio import newickTreeParser
import glob
import os
import subprocess
import sys
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
    print "USAGE: %s --simDir simulationDirectory/ --jobFile JOB_FILE" %( sys.argv[0] )
    print __doc__
    sys.exit(2)

def customOptions(parser):
    parser.add_option('-s', '--simDir',dest='simDir',
                      help='Simulation directory.')
    parser.add_option('-j', '--jobFile',dest='jobFile',
                      help='jobFile, passed in by jobTree.py.')
    parser.add_option('--debug', action='store_true', dest='isDebug',
                      default=False, help='Performs no operations, prints out status of variables.')
    
def customOptionsCheck( options, usage ):
    if options.simDir == None:
        sys.stderr.write('%s: Error, specify --simDir.\n' % sys.argv[0] )
        usage()
    if not os.path.isdir( options.simDir ):
        sys.stderr.write('%s: Error, --simDir "%s" not a directory!\n' % (sys.argv[0], options.simDir))
        usage()
    options.simDir  = os.path.abspath( options.simDir )

    xmlTree = ET.parse( os.path.join( options.simDir, 'simulationInfo.xml' ) )
    tree = xmlTree.find('tree')
    nt = newickTreeParser( tree.text, 0.0 )
    rootNameObj = xmlTree.find('rootDir')
    options.rootName = os.path.basename( rootNameObj.text )
    options.rootDir  = os.path.abspath(  rootNameObj.text )
    if options.jobFile == None:
        sys.stderr.write('%s: Error, specify --jobFile.\n' % sys.argv[0])
        usage()

def cycleDirectoriesOnly( dirs ):
    """directoriesOnly() takes a list of items from a directory
    and creates a new list that only contains desired directories.
    """
    cDirs = []
    for d in dirs:
        if os.path.isdir( d ):
            if os.path.exists( os.path.join( d, 'cycleInfo.xml' ) ):
                cDirs.append( os.path.abspath( d ) )
    return cDirs

class CycleProblem( Exception ):
    pass

def parentCycleDirName( dir ):
    if not os.path.exists( os.path.join( dir, 'cycleInfo.xml' ) ):
        raise CycleProblem( '%s is not a cycleDir!' % dir )
    xmlTree = ET.parse( os.path.join( dir, 'cycleInfo.xml' ) ).getroot()
    if xmlTree is None:
        raise CycleProblem('%s does not have valid xml in cycleInfo.xml' % dir )
    parentDir = xmlTree.find('parentDir')
    if parentDir is None:
        raise CycleProblem('%s does not have tag "parentDir" in cycleInfo.xml' % dir )
    return parentDir.text

def cycleIsBranchOrRoot( dir, options ):
    """ returns True if the parent cycle was either a branch point or the root
    """
    if os.path.samefile( dir, options.rootDir ):
        return True
    file = os.path.join( dir, 'cycleInfo.xml')
    infoTree = ET.parse(file)
    root = infoTree.getroot()
    if len( root.findall('followUpCommand') ) == 3:
        # 3 is a branch, 2 is a stem
        return True
    return False

def issueJobs( dirs, options ):
    numJobs = 0
    jobCMD = ''
    for B in dirs:
        if os.path.basename( B ) == options.rootName:
            continue
        A = parentCycleDirName( B )
        if not os.path.exists( os.path.join( B, 'aln.rev' ) ) and ( os.path.exists( os.path.join( A, 'aln.rev' ) ) or cycleIsBranchOrRoot( A, options ) ):
            jobCMD += CMD_EVAL_BIN +\
                      ' --statXML '+os.path.join(B, 'logs','trans.info.xml')+\
                      ' JOB_FILE "'
            # command_0
            jobCMD +=LSC.commandPacker(TRANS_BIN +\
                                         ' -in1 '+os.path.join(B, 'inter','inter.aln.rev')+ \
                                         ' -in2 '+os.path.join(B, 'intra', 'intra.aln.rev')+ \
                                         ' -out '+os.path.join(B, 'inter-intra.aln.rev')+ \
                                         ' -log '+os.path.join(B, 'logs', 'transalign1.log'))
            # command_1
            if cycleIsBranchOrRoot( A, options ):
                # In these cases the alignment above the branch point should not be carried
                # into the the descendant genomes. Alignments should only go back to the most
                # recent branch point.
                jobCMD += LSC.commandPacker(LINK_BIN +\
                                              ' -s '+os.path.join(B,'inter-intra.aln.rev')+\
                                              ' '+os.path.join(B,'aln.rev'))
            else:
                jobCMD += LSC.commandPacker(TRANS_BIN +\
                                              ' -in1 '+os.path.join(A, 'aln.rev')+ \
                                              ' -in2 '+os.path.join(B, 'inter-intra.aln.rev')+ \
                                              ' -out '+os.path.join(B, 'aln.rev')+ \
                                              ' -log '+os.path.join(B, 'logs', 'transalign2.log'))
            # command_2
            jobCMD += LSC.commandPacker(EVO_BIN+\
                                          ' -cdsalns '+ os.path.join(B, 'intra', 'intra.aln.rev')+ \
                                          ' -alns '   + os.path.join(B, 'stats', 'cds_alns.rev')+ \
                                          ' -annots1 '+ os.path.join(B, 'inter', 'inter.outannots.gff')+ \
                                          ' -annots2 '+ os.path.join(B, 'annots.gff')+ \
                                          ' -log '+ os.path.join(B, 'logs', 'cds_alns.log'))
            # command_3
            jobCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                         ' -fromrev '+os.path.join(B, 'inter', 'inter.aln.rev')+\
                                         ' -tocmap '+os.path.join(B, 'stats', 'img.cycle.cmap.pdf')+\
                                         ' -blocksize '+str(DRAW_REV_BLOCK_SIZE))
            # command_4
            jobCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                         ' -fromrev '+os.path.join(B, 'inter', 'inter.aln.rev')+\
                                         ' -tolmap ' +os.path.join(B, 'stats', 'img.cycle.lmap.png')+\
                                         ' -npp '    +str(DRAW_REV_NT_PER_PIX))
            # command_5
            jobCMD +=LSC.commandPacker(EVO_BIN+\
                                         ' -nologcmdlineandtime '+\
                                         ' -ancstats '+os.path.join(B, 'aln.rev') +\
                                         ' -log '+ os.path.join(B, 'stats', 'tmpstats.root.difflength.tmp'))
            # command_6
            jobCMD +=LSC.commandPacker(EVO_BIN+\
                                         ' -nologcmdlineandtime '+\
                                         ' -ancstats '+os.path.join(B, 'intra', 'intra.aln.rev') +\
                                         ' -log '+ os.path.join(B, 'stats', 'tmpstats.difflength.tmp'))
            # command_7
            jobCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                         ' -fromrev '+os.path.join(B, 'aln.rev') +\
                                         ' -tocmap '+os.path.join(B, 'stats', 'img.root.cmap.pdf') +\
                                         ' -blocksize '+str(DRAW_REV_BLOCK_SIZE))
            # command_8
            jobCMD +=LSC.commandPacker(DRAWREV_BIN+\
                                         ' -fromrev '+os.path.join(B, 'aln.rev')+\
                                         ' -tolmap ' +os.path.join(B, 'stats', 'img.root.lmap.png')+\
                                         ' -npp '    +str(DRAW_REV_NT_PER_PIX))
            # command_9
            jobCMD +=LSC.commandPacker(EVO_BIN+\
                                         ' -cdsalns '+ os.path.join(B, 'aln.rev')+ \
                                         ' -alns '   + os.path.join(B, 'stats', 'cds_alns.rev')+ \
                                         ' -annots1 '+ os.path.join(options.rootDir, 'stats', 'cds_annots.gff')+ \
                                         ' -annots2 '+ os.path.join(B, 'annots.gff')+ \
                                         ' -log '+ os.path.join(B, 'logs', 'root.cds_alns.log'))
            # command_10
            jobCMD +=LSC.commandPacker(EVO_BIN+\
                                         ' -getcodonsubs '+os.path.join(B,'stats', 'cds_alns.rev')+\
                                         ' -out '+os.path.join(B, 'stats', 'codon_sub.txt')+\
                                         ' -log '+os.path.join(B, 'logs','getcodonsubs.log'))
            # command_11
            jobCMD +=LSC.commandPacker(EVO_BIN+\
                                         ' -getcodonsubs '+os.path.join(B,'stats', 'cds_alns.rev')+\
                                         ' -out '+os.path.join(B, 'stats', 'root.codon_sub.txt')+\
                                         ' -log '+os.path.join(B, 'logs','root.getcodonsubs.log'))
            # command_12
            jobCMD +=LSC.commandPacker(CODONREPORT_BIN+\
                                         ' '+os.path.basename(A)+\
                                         ' '+os.path.basename(B)+\
                                         ' < '+os.path.join(B, 'stats', 'codon_sub.txt')+\
                                         ' > '+os.path.join(B, 'stats', 'protstats.txt'))
            # command_13
            jobCMD +=LSC.commandPacker(CODONREPORT_BIN+\
                                         ' '+os.path.basename(options.rootDir)+\
                                         ' '+os.path.basename(B)+\
                                         ' < '+os.path.join(B, 'stats', 'root.codon_sub.txt')+\
                                         ' > '+os.path.join(B, 'stats', 'root.protstats.txt'))
            # command_14
            jobCMD +=LSC.commandPacker(CAT_BIN+\
                                         ' '+ os.path.join(B, 'stats', 'tmpstats.difflength.tmp')+\
                                         ' '+ os.path.join(B, 'stats', 'tmpstats.diffcompost.tmp')+\
                                         ' '+ os.path.join(B,'stats','tmpstats.diffannots.tmp')+\
                                         ' > '+ os.path.join(B, 'stats', 'childParent.diff.txt'))
            # command_15
            jobCMD +=LSC.commandPacker(CAT_BIN+\
                                         ' '+ os.path.join(B, 'stats', 'tmpstats.root.difflength.tmp')+\
                                         ' '+ os.path.join(B, 'stats', 'tmpstats.root.diffcompost.tmp')+\
                                         ' '+ os.path.join(B, 'stats', 'tmpstats.root.diffannots.tmp')+\
                                         ' > '+ os.path.join(B, 'stats', 'childRoot.diff.txt'))
            numJobs += 1
    jobCMD += '"'
    if numJobs:
        xmlTree = ET.parse( options.jobFile )
        jobElm=xmlTree.getroot()
        childrenElm = xmlTree.find('children')
        newChild = ET.SubElement(childrenElm, 'child')
        newChild.attrib['command'] = jobCMD
        followUpCommand = sys.argv[0] +\
                          ' --simDir '+options.simDir+\
                          ' --jobFile JOB_FILE '
        jobElm.attrib['command'] = followUpCommand
        xmlTree.write(options.jobFile)
            

def main():
    parser=OptionParser()
    customOptions(parser)
    (options, args) = parser.parse_args()
    customOptionsCheck( options, usage )

    dirs = glob.glob( os.path.join(options.simDir, '*') )

    dirs = cycleDirectoriesOnly( dirs )

    issueJobs( dirs, options )
    

if __name__ == '__main__':
    main()
    
