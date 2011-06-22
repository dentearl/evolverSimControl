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
import evolverSimControl.lib.libSimControl as lsc
import glob
from jobTree.scriptTree.target import Target
import os
import re
from sonLib.bioio import newickTreeParser
from sonLib.bioio import logger
import subprocess
import sys

class BadInputError(TypeError):
    pass

class SimTree(Target):
    """
    The SimTree class runs the entire simulation. It begins by either calling one Tree() or
    two Tree() targets depending on the newick tree. It then waits for the simulation to complete.
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        nt = newickTreeParser(self.options.inputNewick, 0.0)
        if nt.distance == 0:
            self.addChildTarget(Tree(lsc.tree2str(nt.left), self.options.parentDir,  
                                     'left', self.options))
            self.addChildTarget(Tree(lsc.tree2str(nt.right), self.options.parentDir, 
                                     'right', self.options))
        else:
            self.addChildTarget(Tree(lsc.tree2str(nt), self.options.parentDir, 
                                     'stem', self.options))

class Tree(Target):
    """ The Tree class launches Cycle()'s as children depending on the 
    current tree and issues a follow-on, TreeFollow().
    branchStr is used in TreeFollow to adjust the "random" seed.
    """
    def __init__(self, thisNewickStr, parentDir, branchStr, options):
        Target.__init__(self)
        (self.thisNewickStr, self.thisStepLength)  = lsc.takeNewickStep(thisNewickStr, options)
        self.parentDir = parentDir
        self.thisBranchStr = branchStr # either 'left','right', 'stem'
        self.options = options
        
    def run(self):
        logger.info('Tree object running, %s\n' % self.parentDir)
        if self.thisBranchStr in ['left', 'stem']:
            self.addChildTarget(Transalign(self.parentDir, lsc.getParentDir(self.parentDir), 
                                           self.options))
            self.addChildTarget(Stats(self.parentDir, lsc.getParentDir(self.parentDir), 
                                      self.options))
        self.addChildTarget(Cycle(self.thisNewickStr, self.parentDir, 
                                  self.thisStepLength, self.options))
        self.setFollowOnTarget(TreeFollow(self.thisNewickStr, self.parentDir, 
                                          self.thisBranchStr, self.options))

class TreeFollow(Target):
    """ TreeFollow launches three to four children: Stats and Transalign for the 
    predecessor Tree step and then one or two new Tree steps, depending on whether
    or not the processor was an internal branch point.
    """
    def __init__(self, thisNewickStr, thisGrandParentDir, branchStr, options):
        Target.__init__(self)
        self.thisNewickStr = thisNewickStr
        self.thisGrandParentDir = thisGrandParentDir
        self.options = options
        if self.options.seed != 'stochastic':
            if branchStr == 'left':
                self.options.seed += 47
            elif branchStr == 'right':
                self.options.seed -= 61
            else:
                self.options.seed += 13
            self.options.seed = abs(self.options.seed)

    def run(self):
        logger.info('TreeFollow object running, %s\n' % self.thisGrandParentDir)
        nt = newickTreeParser(self.thisNewickStr, 0.0)
        name = lsc.nameTree(nt)
        commonParentDir = os.path.abspath(os.path.join(self.options.simDir, name))
        if nt.distance == 0:
            if nt.internal:
                # branch point
                branches = { 'left' : lsc.tree2str(nt.left),
                             'right': lsc.tree2str(nt.right) }
                for b in branches:
                    if not lsc.nodeIsLeaf(branches[b]):
                        self.addChildTarget(Tree(branches[b], commonParentDir, b, self.options))
                        childDir = lsc.treeStr2Dir(lsc.takeNewickStep(branches[b], self.options)[0], 
                                                    self.options.simDir)
            else:
                # follow up to leaf cycles... Transalign and Stats only
                self.setFollowOnTarget(LeafCleanUp(commonParentDir, 
                                                   self.thisGrandParentDir, self.options))
        else:
            # stem with distance
            self.addChildTarget(Tree(lsc.tree2str(nt), commonParentDir, 'stem', self.options))
            childDir = lsc.treeStr2Dir(lsc.takeNewickStep(lsc.tree2str(nt), self.options)[0], 
                                        self.options.simDir)

class LeafCleanUp(Target):
    """ LeafCleanUp is called by the TreeFollow object. It only runs
    on leaf cycles and it runs the final Transalign and Stats steps for
    those cycles.
    """
    def __init__(self, thisDir, parentDir, options):
        Target.__init__(self)
        self.thisDir = thisDir
        self.parentDir = parentDir
        self.options = options
    def run(self):
        logger.info('LeafCleanUp object running, %s\n' % self.thisDir)
        self.addChildTarget(Transalign(self.thisDir, self.parentDir, self.options))
        self.addChildTarget(Stats(self.thisDir, self.parentDir, self.options))

class Cycle(Target):
    """ The Cycle class creates the necessary directory structure for the
    given Cycle and then launches CycleStep1 as a child.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepLength, options):
        Target.__init__(self)
        self.thisNewickStr = thisNewickStr
        self.thisParentDir = thisParentDir
        self.thisStepLength  = thisStepLength
        self.options = options
        self.thisDir = lsc.treeStr2Dir(self.thisNewickStr, options.simDir)
        self.theChild  = os.path.basename(self.thisDir)
        self.theParent = os.path.basename(self.thisParentDir)
    def run(self):
        logger.info('Cycle object running, %s\n' % self.thisDir)
        os.mkdir(self.thisDir)
        for d in ['inter', 'intra', 'logs', 'stats', 'xml']:
            if not os.path.exists(os.path.join(self.thisDir, d)):
                os.mkdir(os.path.join(self.thisDir, d))
        if not self.options.noMEs:
            os.mkdir(os.path.join(self.thisDir, 'mobiles'))
        lsc.createNewCycleXmls(self.thisDir, self.thisParentDir, self.thisStepLength, 
                                self.thisNewickStr, self.options)
        self.addChildTarget(CycleStep1(self.thisNewickStr, self.thisParentDir, 
                                       self.thisStepLength, self.options))

class CycleStep1(Cycle):
    """ CycleStep1 consists of an evolver inter step and then the mobiles step.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepLength, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepLength, options)
    def run(self):
        logger.info('CycleStep1 object running, %s\n' % self.thisDir)
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep1_start')
        
        lsc.runEvolverInterCmds(self.thisDir, self.thisParentDir, self.theChild, self.theParent,
                                self.thisStepLength, self.options.seed, self.options.paramsDir,
                                self.getLocalTempDir(), self.options)

        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep1_end')
        self.setFollowOnTarget(CycleStep2(self.thisNewickStr, self.thisParentDir, 
                                            self.thisStepLength, self.options))

class CycleStep2(Cycle):
    """ CycleStep2 sets up the individual evolver intra steps which are run in
    parallel, one per chromosome.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepLength, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepLength, options)
    def run(self):
        logger.info('CycleStep2 object running, %s\n' % self.thisDir)
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep2_start')
        lsc.verifyDirExists(os.path.join(self.thisDir, 'inter'))
        lsc.verifyFileExists(os.path.join(self.thisDir, 'inter', 'inter.chrnames.txt'))
        f = open(os.path.join(self.thisDir, 'inter', 'inter.chrnames.txt'), 'r')
        for chrom in f:
            chrom = chrom.strip()
            self.addChildTarget(CycleStep2Chromosome(self.thisNewickStr, self.thisParentDir,
                                                       self.thisStepLength, chrom, self.options))
        f.close()
        self.setFollowOnTarget(CycleStep3(self.thisNewickStr, self.thisParentDir, 
                                            self.thisStepLength, self.options))

class CycleStep2Chromosome(Cycle):
    """ CycleStep2Chromosome is called by CycleStep2. This corresponds to the 
    evolver intra (within chromosome) step.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepLength, thisChr, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepLength, options)
        self.thisChr = thisChr
    def run(self):
        logger.info('CycleStep2Chromosome object running, %s %s\n' % (self.thisDir, self.thisChr))
        lsc.verifyDirExists(self.thisDir)
        if not os.path.exists(os.path.join(self.thisDir, 'xml', 'cycle.%s.xml' % self.thisChr)):
            lsc.newInfoXml(os.path.join(self.thisDir, 'xml', 'cycle.%s.xml' % self.thisChr))
            lsc.addTimestampsTag(os.path.join(self.thisDir, 'xml', 'cycle.%s.xml' % self.thisChr))
            lsc.subTypeTimestamp(self.thisDir, 'cycleChr', 
                                 'CycleStep2Chr_%s_start' % self.thisChr, self.thisChr)

        # evolver intra on one chromosome
        cmds = lsc.evolverIntraStepCmd(self.thisDir, self.theChild, self.thisStepLength, 
                                       self.thisChr, self.options.seed, 
                                       self.options.paramsDir, self.getLocalTempDir(), self.options)
        lsc.runCommands(cmds, self.getLocalTempDir())

        # evolver conversion from .rev to fasta in localTempDir
        cmds = lsc.evolverIntraStepToFastaCmd(self.thisDir, self.thisStepLength, self.thisChr, 
                                              self.options.paramsDir, self.getLocalTempDir())
        lsc.runCommands(cmds, self.getLocalTempDir())
            
        # trf wrapper
        lsc.callEvolverIntraStepTRFCmd(self.thisDir, self.thisChr, self.getLocalTempDir())
        
        # move the resulting trf files out of localTempDir
        cmds = lsc.evolverIntraStepMoveTRFCmd(self.thisDir, self.thisChr, self.getLocalTempDir())
        lsc.runCommands(cmds, self.getLocalTempDir(), mode='p')
        
        lsc.subTypeTimestamp(self.thisDir, 'cycleChr', 
                             'CycleStep2Chr_%s_end' % self.thisChr, self.thisChr)
        lsc.addEndTimeAttribute(os.path.join(self.thisDir, 'xml', 'cycle.%s.xml' % self.thisChr))

class CycleStep3(Cycle):
    """ CycleStep3 
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepLength, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepLength, options)
    def run(self):
        logger.info('CycleStep3 object running, %s\n' % self.thisDir)
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep2_end')
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep3_start')

        # trfBig
        # lsc.runMergeTrfBedsToGff(self.thisDir)

        # trf
        regex = r'^(chr\S+)\.outseq\.fa.*\.dat'
        pat = re.compile(regex)
        files = glob.glob(os.path.join(self.thisDir, 'intra', '*.dat'))
        cmds = []
        outPipes = []
        followCmds = []
        followPipes = []
        for f in files:
            # each file is the trf output for one chromosome
            m = re.match(regex, os.path.basename(f))
            if m is None:
                raise RuntimeError('Regex "%s" failed on filename %s' % (regex, os.path.basename(f)))
            outname = os.path.join(self.thisDir, 'intra', m.group(1) + 'trfannots.gff')
            if not os.path.exists(outname):
                # convert the .dat to .gff
                cmd = [lsc.which('evolver_trf2gff.py'), f]
                cmds.append(cmd)
                outPipes.append(outname + '.tmp')
                # atomic files
                followCmds.append([lsc.which('mv'), outname + '.tmp', outname])
                followPipes.append(None)
        lsc.runCommands(cmds, self.getLocalTempDir(), outPipes = outPipes, mode='p')
        lsc.runCommands(followCmds, self.getLocalTempDir(), outPipes = followPipes, mode='p')
        
        catCmd, evoCmd, cvtCmd, followCmds = lsc.evolverIntraMergeCmds(self.thisDir, self.theChild)
        
        lsc.runCommands([catCmd, evoCmd, cvtCmd], self.getLocalTempDir(),
                         outPipes = [os.path.join(self.thisDir, 'intra', 'evannots.gff.tmp'), None, None], 
                         mode='p')
        lsc.runCommands(followCmds, self.getLocalTempDir())
                
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep3_end')
        self.setFollowOnTarget(CycleStep4(self.thisNewickStr, self.thisParentDir,
                                            self.thisStepLength, self.options))

class CycleStep4(Cycle):
    """ CycleStep4 
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepLength, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepLength, options)
    def run(self):
        logger.info('CycleStep4 object running, %s\n' % self.thisDir)
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep4_start')

        outname = os.path.join(self.thisDir, 'logs', 'gene_deactivation.log')
        if not os.path.exists(outname):
            if not self.options.noGeneDeactivation:
                # by default gene deactivation is turned on.
                cmd = lsc.evolverGeneDeactivationStep(self.thisDir, self.thisParentDir)
                p = subprocess.Popen(cmd, cwd=self.getLocalTempDir(), 
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                out = p.communicate()[0]
                f=open(outname + '.tmp', 'w')
                f.write(out)
                f.close()
                os.rename(outname + '.tmp', outname)
            else:
                # this could cause a proliferation of gene creation.
                cmd = [lsc.which('cp')]
                cmd.append(os.path.join(thisDir, 'intra', 'evannots.gff'))
                cmd.append(os.path.join(thisDir, 'annots.gff'))
                cmds = [cmd]
                cmds.append([lsc.which('touch'), outname])
                lsc.runCommands(cmds, self.getLocalTempDir())
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep4_end')
        lsc.typeTimestamp(os.path.join(self.thisDir), 'cycle', 'end')
        lsc.addEndTimeAttribute(os.path.join(self.thisDir, 'xml', 'cycle.xml'))

class Stats(Target):
    """ The Stats object is a convenience class that launches
    StatsStep1 as a child.
    """
    def __init__(self, thisDir, thisParentDir, options):
        Target.__init__(self)
        self.thisDir = thisDir
        self.thisParentDir = thisParentDir
        self.options = options
    def run(self):
        if self.thisParentDir is None:
            # happens when thisParentDir is the root
            return
        lsc.verifyDirExists(self.thisDir)
        lsc.newInfoXml(os.path.join(self.thisDir, 'xml', 'stats.xml'))
        lsc.typeTimestamp(os.path.join(self.thisDir), 'stats', 'start')


        self.addChildTarget(StatsStep1(self.thisDir, self.thisParentDir, self.options))

class StatsStep1(Stats):
    """ StatsStep1 
    """
    def __init__(self, thisDir, thisParentDir, options):
        Stats.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        logger.info('StatsStep1 object running, %s\n' % self.thisDir)
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep1_start')
        
        cmds, followCmds, outPipes = lsc.statsStep1CmdsP(self.thisDir, self.thisParentDir)
        lsc.runCommands(cmds, self.getLocalTempDir(), outPipes = outPipes, mode='p')
        lsc.runCommands(followCmds, self.getLocalTempDir())
        cmds, outPipes = lsc.statsStep1CmdsS(self.thisDir, self.thisParentDir)
        lsc.runCommands(cmds, self.getLocalTempDir(), outPipes = outPipes)
        
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep1_end')
        self.setFollowOnTarget(StatsStep2(self.thisDir, self.thisParentDir, self.options))

class StatsStep2(Stats):
    """ StatsStep2
    """
    def __init__(self, thisDir, thisParentDir, options):
        Stats.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        logger.info('StatsStep2 object running, %s\n' % self.thisDir)
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep2_start')

        cmds, outPipes = lsc.statsStep2Cmds(self.thisDir, self.thisParentDir, self.options)
        lsc.runCommands(cmds, self.getLocalTempDir(), outPipes = outPipes)

        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep2_end')
        self.setFollowOnTarget(StatsStep3(self.thisDir, self.thisParentDir, self.options))

class StatsStep3(Stats):
    """ StatsStep3 
    """
    def __init__(self, thisDir, thisParentDir, options):
        Stats.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        logger.info('StatsStep3 object running, %s\n' % self.thisDir)
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep3_start')

        cmds, pipes = lsc.statsStep3Cmds(self.thisDir, self.thisParentDir, self.options)
        lsc.runCommands(cmds, self.getLocalTempDir(), outPipes = pipes)

        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep3_end')
        self.setFollowOnTarget(StatsStep4(self.thisDir, self.thisParentDir, self.options))

class StatsStep4(Stats):
    """ StatsStep4 
    """
    def __init__(self, thisDir, thisParentDir, options):
        Stats.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        logger.info('StatsStep4 object running, %s\n' % self.thisDir)
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep4_start')
        
        cmds, pipes = lsc.statsStep4Cmds(self.thisDir, self.thisParentDir, self.options)
        lsc.runCommands(cmds, self.getLocalTempDir(), outPipes = pipes)

        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep4_end')
        lsc.typeTimestamp(os.path.join(self.thisDir), 'stats', 'end')
        lsc.addEndTimeAttribute(os.path.join(self.thisDir, 'xml', 'stats.xml'))
        lsc.lastOneOutTurnOffTheLightsCycle(self.thisDir)
        if lsc.isLeaf(self.thisDir):
            lsc.lastOneOutTurnOffTheLightsSimulation(self.options.simDir, self.options)

class Transalign(Target):
    """ The Transalign class is a convenience class that
    launches the two TransalignStep classes.
    """
    def __init__(self, thisDir, thisParentDir, options):
        Target.__init__(self)
        self.thisDir = thisDir
        self.thisParentDir = thisParentDir
        self.options = options
    def run(self):
        logger.info('Transalign object running, thisDir: %s thisParentDir: %s\n' 
                    % (self.thisDir, self.thisParentDir))
        if self.thisParentDir is None:
            # happens when thisParentDir is the root
            return
        lsc.verifyDirExists(self.thisDir)
        lsc.newInfoXml(os.path.join(self.thisDir, 'xml', 'transalign.xml'))
        lsc.typeTimestamp(self.thisDir, 'transalign', 'start')

        self.addChildTarget(TransalignStep1(self.thisDir, self.thisParentDir, self.options))

class TransalignStep1(Transalign):
    """ TransalignStep1
    """
    def __init__(self, thisDir, thisParentDir, options):
        Transalign.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        logger.info('TransalignStep1 object running, thisDir: %s thisParentDir: %s\n' 
                    % (self.thisDir, self.thisParentDir))
        lsc.verifyDirExists(self.thisDir)
        lsc.subTypeTimestamp(self.thisDir, 'transalign', 'TransalignStep1_start')
        
        cmds, pipes = lsc.transalignStep1Cmds_1(self.thisDir, self.thisParentDir, self.options)
        lsc.runCommands(cmds, self.getLocalTempDir(), outPipes = pipes)
        
        lsc.runTransalignStep1Cmds_2(self.thisDir, self.thisParentDir, 
                                     self.getLocalTempDir(), self.options)
        
        lsc.subTypeTimestamp(self.thisDir, 'transalign', 'TransalignStep1_end')
        lsc.typeTimestamp(os.path.join(self.thisDir), 'transalign', 'end')
        lsc.addEndTimeAttribute(os.path.join(self.thisDir, 'xml', 'transalign.xml'))
        lsc.lastOneOutTurnOffTheLightsCycle(self.thisDir)
        if lsc.isLeaf(self.thisDir):
            lsc.lastOneOutTurnOffTheLightsSimulation(self.options.simDir, self.options)

class Node:
    """Nodes have one parent and two children,
    unless they are children in which case their
    children list is empty, or they are the root
    node in which case their parent is None.
    This class is used by simCtrl_postSimMafExtractor.py
    """
    def __init__(self):
        self.parent   = None
        self.name     = ''
        self.children = []
        self.isLeaf   = False

class MasterMafGenerator(Target):
    """
    The MasterMafGenerator class runs the entire extraction and merge process. It begins
    by calling the extractionManager which extracts all necessary pairwise mafs.
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        self.setFollowOnTarget(ExtractionManager(self.options))

class ExtractionManager(Target):
    """
    The ExtractionManager class runs the extraction process. 
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        logger.info('ExtractionManager object running, rootDir: %s\n' % (self.options.rootDir))
        nt = newickTreeParser(self.options.inputNewick, 0.0)
        nodesList = []
        leafsDict = {}
        lsc.extractLeafsFromNewick(nt, leafsDict)
        nt.iD = os.path.basename(self.options.rootDir)
        lsc.buildNodesListFromNewick(nt, nodesList, leafsDict)
        if (os.path.exists(os.path.join(self.options.rootDir, 'aln.rev')) and not
            os.path.exists(os.path.join(self.options.rootDir, 'burnin.tmp.maf'))):
            self.addChildTarget(Extract(self.options.rootDir, 'burnin', False, self.options))
        for n in nodesList:
            # parent nodes
            for c in n.children:
                # the child alignment is named for the parent node
                self.addChildTarget(Extract(os.path.join(self.options.simDir, c), n.name,
                                            c in leafsDict, self.options))
        self.setFollowOnTarget(MergeManager(nodesList, leafsDict, self.options))

class Extract(Target):
    """
    The Extract class runs a single extraction. 
    """
    def __init__(self, thisDir, alignName, isLeaf, options):
        Target.__init__(self)
        self.thisDir = thisDir
        self.alignName = alignName
        self.isLeaf = isLeaf
        self.options = options

    def run(self):
        logger.info('Extract object running, thisDir: %s\n' % (self.thisDir))
        if self.isLeaf:
            ext = '.maf'
        else:
            ext = '.tmp.maf'
        outname = os.path.join(self.thisDir, self.alignName + ext)
        if not os.path.exists(outname):
            cmd = [lsc.which('evolver_cvt')]
            cmd.append('-fromrev')
            cmd.append(os.path.join(self.thisDir, 'aln.rev'))
            cmd.append('-tomaf')
            cmd.append(outname + '.tmp')
            cmds = [cmd]
            cmds.append([lsc.which('mv'), outname + '.tmp', outname])
            lsc.runCommands(cmds, self.getLocalTempDir())

class MergeManager(Target):
    """
    The MergeManager class runs the merge process. 
    
    MergeManager takes the nodesList and goes through the
    steps of progressively merging the MAFs, starting from the leaves
    and working its way up the tree.
    ((a, b)E,(c,d)F)root;
    which was earlier decomposed into the:
    aE.maf
    bE.maf
    cF.maf
    dF.maf
    Eroot.maf
    Froot.maf
    would now be combined in the following pattern:
    aE.maf + bE.maf -> abE.maf
    cF.maf + dF.maf -> cdF.maf
    abE.maf + Eroot.maf -> abEroot.maf
    cdF.maf + Froot.maf -> cdFroot.maf
    abEroot.maf + cdFroot.maf -> abcdEFroot.maf
    ... and then we drink root beers and have highfives.
    OR, in the specific directory based data structure
    that we actually use,
    a/E.maf
    b/E.maf
    c/F.maf
    d/F.maf
    E/root.tmp.maf
    F/root.tmp.maf
    a/E.maf b/E.maf -> E/E.maf
    E/E.maf E/root.tmp.maf -> E/root.maf
    F/F.maf F/root.tmp.maf -> F/root.maf
    E/root.maf F/root.maf -> root/all.maf
    OR, most generally,
    child0/parent.maf child1/parent.maf -> parent/parent.maf
    parent/parent.maf parent/grandParent.tmp.maf -> parent/grandParent.maf
    
    """
    def __init__(self, nodesList, leafsDict, options):
        Target.__init__(self)
        self.options = options
        self.nodesList = nodesList
        self.nodeParentDict = lsc.buildNodeParentDict(self.nodesList)
        self.leafsDict = leafsDict
        self.nodeDict = lsc.buildNodesDict(self.nodesList, self.leafsDict)

    def run(self):
        logger.info('Extract object running, rootDir: %s\n' % (self.options.rootDir))
        nt = newickTreeParser(self.options.inputNewick, 0.0)
        nt.iD = os.path.basename(self.options.rootDir)
        self.addChildTarget(MergeTree(nt, self.nodeDict, self.nodeParentDict, 
                                      self.leafsDict, self.options))
        if not self.options.noBurninMerge:
            self.setFollowOnTarget(MergeTreeFollow(nt, self.nodeDict, self.nodeParentDict,
                                                   self.leafsDict, self.options))

class MergeTree(Target):
    """
    The MergeTrees class creates a tree structure of maf merges based on the phylogenetic 
    tree of the simulation.
    
    """
    def __init__(self, nt, nodeDict, nodeParentDict, leafsDict, options):
        Target.__init__(self)
        nt.distance = 0.0
        self.nt = nt
        self.name = lsc.nameTree(self.nt)
        self.nodeDict = nodeDict
        self.nodeParentDict = nodeParentDict
        self.leafsDict = leafsDict
        self.options = options
        if self.name != self.options.rootName:
            self.nodeParent = self.nodeParentDict[self.name]
        else:
            self.nodeParent = self.options.rootName

    def run(self):
        logger.info('MergeTree object running, name: %s\n' % (self.name))
        if self.nt is None:
            return
        for t in [self.nt.left, self.nt.right]:
            if t is None:
                continue
            t.distance = 0.0
            if not os.path.exists(os.path.join(self.options.simDir, lsc.nameTree(t), self.name + '.maf')):
                sys.stderr.write('%s does not exist, recurse\n' 
                                 % os.path.join(self.options.simDir, lsc.nameTree(t), self.name+ '.maf'))
                self.addChildTarget(MergeTree(t, self.nodeDict, 
                                              self.nodeParentDict, self.leafsDict, self.options))
            else:
                sys.stderr.write('%s does exist, don\'t recurse\n' 
                                 % os.path.join(self.options.simDir, lsc.nameTree(t),self.name + '.maf'))
        self.setFollowOnTarget(MergeMafsDown(self.nt, self.nodeDict, self.nodeParentDict, 
                                             self.leafsDict, self.nodeParent, self.options))

class MergeTreeFollow(Target):
    """
    The MergeTreeFollow class checks the rootDir for evidence of a burnin and if such evidence exists,
    it then performs a final merge using that burnin and the current full simulation maf.
    """
    def __init__(self, nt, nodeDict, nodeParentDict, leafsDict, options):
        Target.__init__(self)
        nt.distance = 0.0
        self.nt = nt
        self.name = lsc.nameTree(self.nt)
        self.nodeDict = nodeDict
        self.nodeParentDict = nodeParentDict
        self.leafsDict = leafsDict
        self.options = options
        if self.name != self.options.rootName:
            self.nodeParent = self.nodeParentDict[self.name]
        else:
            self.nodeParent = self.options.rootName
    def run(self):
        logger.info('MergeTreeFollow object running, name: %s\n' % (self.name))
        outname = os.path.join(self.options.rootDir, 'burnin.maf')
        if (os.path.exists(os.path.join(self.options.rootDir, 'burnin.tmp.maf')) 
            and not os.path.exists(outname)):
            treelessRootCmd = ['-treelessRoot2=%s' % lsc.burninRootName(self.options)]
            maf1 = os.path.join(self.options.rootDir, self.options.rootName + '.maf')
            maf2 = os.path.join(self.options.rootDir, 'burnin.tmp.maf')
            drop = os.path.join(self.options.rootDir, 'burnin.dropped.maf')
            cmds = lsc.buildMergeCommand(maf1, maf2, outname, treelessRootCmd, 
                                         self.name, self.options, drop)
            lsc.runCommands(cmds, self.getLocalTempDir())

class MergeMafsDown(MergeTree):
    """
    The MergeMafsDown class runs the merging of two children mafs plus the current maf.
    """
    def __init__(self, nt, nodeDict, nodeParentDict, leafsDict, nodeParent, options):
        Target.__init__(self)
        MergeTree.__init__(self, nt, nodeDict, nodeParentDict, leafsDict, options)
        self.nodeParent = nodeParent

    def run(self):
        logger.info('MergeTreeDown object running, name: %s nodeParent: %s\n' 
                    % (self.name, self.nodeParent))
        treelessRootCmd = []
        for i in xrange(0,2):
            if self.nodeDict[self.name].children[i] in self.leafsDict:
                treelessRootCmd.append('-treelessRoot%d=%s' % (i + 1, self.name))
        ##############################
        # the 'lookdown' aspect of the merge is performed for every node, including the root.
        outname = os.path.join(self.options.simDir, self.name, self.name + '.maf')
        if not os.path.exists(outname):
            maf1 = os.path.join(self.options.simDir, self.nodeDict[self.name].children[0], 
                                self.name + '.maf')
            maf2 = os.path.join(self.options.simDir, self.nodeDict[self.name].children[1], 
                                self.name + '.maf')
            drop = os.path.join(self.options.simDir, self.name, self.name + '.dropped.tab')
            cmds = lsc.buildMergeCommand(maf1, maf2, outname, treelessRootCmd, self.name, 
                                         self.options, drop)
            lsc.runCommands(cmds, self.getLocalTempDir())
        self.setFollowOnTarget(MergeMafsUp(self.nt, self.nodeDict, self.nodeParentDict, 
                                           self.leafsDict, self.nodeParent, self.options))

class MergeMafsUp(MergeTree):
    """
    The MergeMafsUp class runs the second part of a maf merge, merging a maf containing 
    the entire tree including thisDir into the parent of thisDir.
    """
    def __init__(self, nt, nodeDict, nodeParentDict, leafsDict, nodeParent, options):
        Target.__init__(self)
        MergeTree.__init__(self, nt, nodeDict, nodeParentDict, leafsDict, options)
        self.nodeParent = nodeParent

    def run(self):
        logger.info('MergeMafsUp object running, name: %s nodeParent: %s\n' 
                    % (self.name, self.nodeParent))
        ##############################
        # The 'lookup' aspect of the merge is only performed when we are not at the root
        # This merge merges the results of the 'lookdown' merge, that is to say the maf that contains
        # all descendant sequences including the node, with the node-parent maf, to produce a maf
        # that the parent can use to merge its children.
        if self.name == self.options.rootName:
            return
        outname = os.path.join(self.options.simDir, self.name, self.nodeParent + '.maf')
        if not os.path.exists(outname):
            treelessRootCmd = ['-treelessRoot2=%s' % self.nodeParent]
            maf1 = os.path.join(self.options.simDir, self.name, self.name + '.maf')
            maf2 = os.path.join(self.options.simDir, self.name, self.nodeParent + '.tmp.maf')
            drop = os.path.join(self.options.simDir, self.name, self.nodeParent + '.dropped.tab')
            cmds = lsc.buildMergeCommand(maf1, maf2, outname, treelessRootCmd, 
                                         self.name, self.options, drop)
            lsc.runCommands(cmds, self.getLocalTempDir())
    
