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
from jobTree.scriptTree.target import Target
import evolverSimControl.lib.libSimControl as lsc
import os
from sonLib.bioio import newickTreeParser
import subprocess
import sys

class BadInputError(ValueError):
    pass
class ProgramDoesNotExistError(ValueError):
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
        (self.thisNewickStr, self.thisStepSize)  = lsc.takeNewickStep(thisNewickStr, options)
        self.parentDir = parentDir
        self.thisBranchStr = branchStr # either 'left','right', 'stem'
        self.options = options
        
    def run(self):
        if self.thisBranchStr in ['left', 'stem']:
            self.addChildTarget(Transalign(self.parentDir, lsc.getParentDir(self.parentDir), 
                                           self.options))
            self.addChildTarget(Stats(self.parentDir, lsc.getParentDir(self.parentDir), 
                                      self.options))
        self.addChildTarget(Cycle(self.thisNewickStr, self.parentDir, 
                                  self.thisStepSize, self.options))
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
                self.setFollowOnTarget(LeafCleanUp(commonParentDir, self.thisGrandParentDir, self.options))
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
        self.addChildTarget(Transalign(self.thisDir, self.parentDir, self.options))
        self.addChildTarget(Stats(self.thisDir, self.parentDir, self.options))

class Cycle(Target):
    """ The Cycle class creates the necessary directory structure for the
    given Cycle and then launches CycleStep1 as a child.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Target.__init__(self)
        self.thisNewickStr = thisNewickStr
        self.thisParentDir = thisParentDir
        self.thisStepSize  = thisStepSize
        self.options = options
        self.thisDir = lsc.treeStr2Dir(self.thisNewickStr, options.simDir)
        self.theChild  = os.path.basename(self.thisDir)
        self.theParent = os.path.basename(self.thisParentDir)
    def run(self):
        os.mkdir(self.thisDir)
        for d in ['inter', 'intra', 'logs', 'mobiles', 'stats', 'xml']:
            os.mkdir(os.path.join(self.thisDir, d))
        lsc.createNewCycleXmls(self.thisDir, self.thisParentDir, self.thisStepSize, 
                                self.thisNewickStr, self.options)
        self.addChildTarget(CycleStep1(self.thisNewickStr, self.thisParentDir, 
                                       self.thisStepSize, self.options))

class CycleStep1(Cycle):
    """ CycleStep1 consists of an evolver inter step and then the mobiles step.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep1_start')
        lsc.verifyDirExists(self.thisDir)
        cmds = []
        cmd = lsc.evolverInterStepCmd(self.thisDir, self.thisParentDir, self.theChild, 
                                       self.thisStepSize, self.options.seed, self.options.paramsDir)
        cmds.append(cmd)
        p1 = subprocess.Popen(cmd, cwd=self.getLocalTempDir())
        if not self.options.noMEs:
            cmd = lsc.evolverInterStepMobilesCmd(self.thisDir, self.thisParentDir, self.theParent, 
                                                  self.thisStepSize, self.options.paramsDir)
            cmds.append(cmd)
            p2 = subprocess.Popen(cmd, cwd=self.getLocalTempDir(), stdout=subprocess.PIPE)
            f = open(os.path.join(self.thisDir, 'logs', 'mobiles.log'), 'w')
            f.write(p2.communicate()[0])
            f.close()
            lsc.handleReturnCode(p2.returncode, cmds[1])
            mvCmds= lsc.evolverInterStepMobilesMoveCmd(self.getLocalTempDir(), self.thisDir)
            lsc.runCommands(mvCmds, self.getLocalTempDir(), mode='p')
        p1.wait()
        lsc.handleReturnCode(p1.returncode, cmds[0])
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep1_end')
        self.setFollowOnTarget(CycleStep2(self.thisNewickStr, self.thisParentDir, 
                                            self.thisStepSize, self.options))

class CycleStep2(Cycle):
    """ CycleStep2 sets up the individual evolver intra steps which are run in
    parallel, one per chromosome.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep2_start')
        lsc.verifyDirExists(self.thisDir)
        lsc.verifyDirExists(os.path.join(self.thisDir, 'inter'))
        lsc.verifyFileExists(os.path.join(self.thisDir, 'inter', 'inter.chrnames.txt'))
        f = open(os.path.join(self.thisDir, 'inter', 'inter.chrnames.txt'), 'r')
        for chrom in f:
            chrom = chrom.strip()
            self.addChildTarget(CycleStep2Chromosome(self.thisNewickStr, self.thisParentDir,
                                                       self.thisStepSize, chrom, self.options))
        f.close()
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep2_end')
        
        self.setFollowOnTarget(CycleStep3(self.thisNewickStr, self.thisParentDir, 
                                            self.thisStepSize, self.options))

class CycleStep2Chromosome(Cycle):
    """ CycleStep2Chromosome is called by CycleStep2. This corresponds to the 
    evolver intra (within chromosome) step.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, thisChr, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
        self.thisChr = thisChr
    def run(self):
        lsc.newInfoXml(os.path.join(self.thisDir, 'xml', 'cycle.%s.xml' % self.thisChr))
        lsc.addTimestampsTag(os.path.join(self.thisDir, 'xml', 'cycle.%s.xml' % self.thisChr))
        lsc.subTypeTimestamp(self.thisDir, 'cycleChr', 'CycleStep2Chr_%s_start' % self.thisChr, self.thisChr)
        lsc.verifyDirExists(self.thisDir)
        cmds = []
        # evolver intra on one chromosome
        cmd = lsc.evolverIntraStepCmd(self.thisDir, self.theChild, self.thisStepSize, 
                                       self.thisChr, self.options.seed, 
                                       self.options.paramsDir, self.getLocalTempDir(), self.options)
        cmds.append(cmd)
        lsc.runCommands(cmds, self.getLocalTempDir())
        cmds = []
        # evolver conversion from .rev to fasta in localTempDir
        cmd = lsc.evolverIntraStepToFastaCmd(self.thisDir, self.thisStepSize, self.thisChr, 
                                              self.options.paramsDir, self.getLocalTempDir())
        cmds.append(cmd)
        lsc.runCommands(cmds, self.getLocalTempDir())
            
        # trf wrapper
        lsc.callEvolverIntraStepTRFCmd(self.thisDir, self.thisChr, self.getLocalTempDir())
        
        # move the resulting trf files out of localTempDir
        cmds = lsc.evolverIntraStepMoveTRFCmd(self.thisDir, self.thisChr, self.getLocalTempDir())
        lsc.runCommands(cmds, self.getLocalTempDir(), mode='p')
        
        lsc.subTypeTimestamp(self.thisDir, 'cycleChr', 'CycleStep2Chr_%s_end' % self.thisChr, self.thisChr)

class CycleStep3(Cycle):
    """ CycleStep3 
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep3_start')
        lsc.verifyDirExists(self.thisDir)
        
        # cmd = [lsc.which('evolver_trf2gff.py')]
        # files = glob.glob(os.path.join(self.thisDir, 'intra', '*.dat'))
        # for f in files:
        #     cmd.append(f)
        # lsc.runCommands([cmd], self.getLocalTempDir(), [os.path.join(self.thisDir, 'intra', 'trfannots.gff')])
        lsc.runMergeTrfBedsToGff(self.thisDir)
        
        catCmd, evoCmd, cvtCmd = lsc.evolverIntraMergeCmds(self.thisDir, self.theChild)
        
        lsc.runCommands([catCmd, evoCmd, cvtCmd], self.getLocalTempDir(),
                         [os.path.join(self.thisDir, 'intra', 'evannots.gff'), None, None], 
                         mode='p')
                
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep3_end')
        self.setFollowOnTarget(CycleStep4(self.thisNewickStr, self.thisParentDir,
                                            self.thisStepSize, self.options))

class CycleStep4(Cycle):
    """ CycleStep4 
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep4_start')
        lsc.verifyDirExists(self.thisDir)
        if not self.options.noGeneDeactivation:
            # by default gene deactivation is turned on.
            cmd = lsc.evolverGeneDeactivationStep(self.thisDir, self.thisParentDir)
            p = subprocess.Popen(cmd, cwd=self.getLocalTempDir(), 
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            out = p.communicate()[0]
            f=open(os.path.join(self.thisDir, 'logs', 'gene_deactivation.log'), 'w')
            f.write(out)
            f.close()
        else:
            # this could cause a proliferation of gene creation.
            cmd = [lsc.which('cp')]
            cmd.append(os.path.join(thisDir, 'intra', 'evannots.gff'))
            cmd.append(os.path.join(thisDir, 'annots.gff'))
            p = subprocess.Popen(cmd, cwd=self.getLocalTempDir())
            lsc.handleReturnCode(p.returncode, cmd)
        lsc.subTypeTimestamp(self.thisDir, 'cycle', 'CycleStep4_end')
        lsc.typeTimestamp(os.path.join(self.thisDir), 'cycle', 'end')

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
        lsc.newInfoXml(os.path.join(self.thisDir, 'xml', 'stats.xml'))
        lsc.typeTimestamp(os.path.join(self.thisDir), 'stats', 'start')
        lsc.verifyDirExists(self.thisDir)

        self.addChildTarget(StatsStep1(self.thisDir, self.thisParentDir, self.options))

class StatsStep1(Stats):
    """ StatsStep1 
    """
    def __init__(self, thisDir, thisParentDir, options):
        Stats.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep1_start')
        lsc.verifyDirExists(self.thisDir)
        
        cmds, pipes = lsc.statsStep1CmdsP(self.thisDir, self.thisParentDir)
        lsc.runCommands(cmds, self.getLocalTempDir(), pipes, mode='p')
        cmds, pipes = lsc.statsStep1CmdsS(self.thisDir, self.thisParentDir)
        lsc.runCommands(cmds, self.getLocalTempDir(), pipes)
        
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep1_end')
        self.setFollowOnTarget(StatsStep2(self.thisDir, self.thisParentDir, self.options))

class StatsStep2(Stats):
    """ StatsStep2
    """
    def __init__(self, thisDir, thisParentDir, options):
        Stats.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep2_start')
        lsc.verifyDirExists(self.thisDir)

        cmds, pipes = lsc.statsStep2Cmds(self.thisDir, self.thisParentDir, self.options)
        lsc.runCommands(cmds, self.getLocalTempDir(), pipes)

        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep2_end')
        self.setFollowOnTarget(StatsStep3(self.thisDir, self.thisParentDir, self.options))

class StatsStep3(Stats):
    """ StatsStep3 
    """
    def __init__(self, thisDir, thisParentDir, options):
        Stats.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep3_start')
        lsc.verifyDirExists(self.thisDir)

        cmds, pipes = lsc.statsStep3Cmds(self.thisDir, self.thisParentDir, self.options)
        lsc.runCommands(cmds, self.getLocalTempDir(), pipes)

        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep3_end')
        self.setFollowOnTarget(StatsStep4(self.thisDir, self.thisParentDir, self.options))

class StatsStep4(Stats):
    """ StatsStep4 
    """
    def __init__(self, thisDir, thisParentDir, options):
        Stats.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'stats', 'StatsStep4_start')
        lsc.verifyDirExists(self.thisDir)
        
        cmds, pipes = lsc.statsStep4Cmds(self.thisDir, self.thisParentDir, self.options)
        lsc.runCommands(cmds, self.getLocalTempDir(), pipes)

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
        if self.thisParentDir is None:
            # happens when thisParentDir is the root
            return
        lsc.newInfoXml(os.path.join(self.thisDir, 'xml', 'transalign.xml'))
        lsc.typeTimestamp(self.thisDir, 'transalign', 'start')
        lsc.verifyDirExists(self.thisDir)

        self.addChildTarget(TransalignStep1(self.thisDir, self.thisParentDir, self.options))

class TransalignStep1(Transalign):
    """ TransalignStep1
    """
    def __init__(self, thisDir, thisParentDir, options):
        Transalign.__init__(self, thisDir, thisParentDir, options)
    def run(self):
        lsc.subTypeTimestamp(self.thisDir, 'transalign', 'TransalignStep1_start')
        lsc.verifyDirExists(self.thisDir)
        
        cmds, pipes = lsc.transalignStep1Cmds_1(self.thisDir, self.thisParentDir, self.options)
        lsc.runCommands(cmds, self.getLocalTempDir(), pipes)
        
        lsc.runTransalignStep1Cmds_2(self.thisDir, self.thisParentDir, self.getLocalTempDir(), self.options)
        
        lsc.subTypeTimestamp(self.thisDir, 'transalign', 'TransalignStep1_end')
        lsc.typeTimestamp(os.path.join(self.thisDir), 'transalign', 'end')
        lsc.addEndTimeAttribute(os.path.join(self.thisDir, 'xml', 'transalign.xml'))
        lsc.lastOneOutTurnOffTheLightsCycle(self.thisDir)
        if lsc.isLeaf(self.thisDir):
            lsc.lastOneOutTurnOffTheLightsSimulation(self.options.simDir, self.options)
