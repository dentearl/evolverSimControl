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
import sys

class BadInputError( ValueError ):
    pass
class ProgramDoesNotExistError( ValueError ):
    pass

class SimTree( Target ):
    """
    The SimTree class runs the entire simulation. It begins by either calling one Tree() or
    two Tree() targets depending on the newick tree. It then waits for the simulation to complete.
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        nt = newickTreeParser( self.options.inputNewick, 0.0 )
        if nt.distance == 0:
            self.addChildTarget( Tree( lsc.tree2str(nt.left), self.options.parentDir,  'left', self.options) )
            self.addChildTarget( Tree( lsc.tree2str(nt.right), self.options.parentDir, 'right', self.options) )
        else:
            self.addChildTarget( Tree( lsc.tree2str(nt), self.options.parentDir, 'stem', self.options) )

class Tree( Target ):
    """ The Tree class launches Cycle()'s as children depending on the 
    current tree and issues a follow-on, TreeFollow().
    branchStr is used in TreeFollow to adjust the "random" seed.
    """
    def __init__(self, thisNewickStr, parentDir, branchStr, options):
        Target.__init__(self)
        (self.thisNewickStr, self.thisStepSize)  = lsc.takeNewickStep( thisNewickStr, options )
        self.parentDir = parentDir
        self.thisBranchStr = branchStr # either 'left','right', 'stem'
        self.options = options
        
    def run(self):
        lsc.myLog('Tree.run(), %s %s\n' % (self.thisNewickStr, self.thisStepSize))
        self.addChildTarget( Cycle( self.thisNewickStr, self.parentDir, self.thisStepSize, self.options ))
        self.setFollowOnTarget( TreeFollow( self.thisNewickStr, self.thisBranchStr, self.options))

class TreeFollow( Target ):
    """ TreeFollow launches three to four children: Stats and Transalign for the 
    predecessor Tree step and then one or two new Tree steps, depending on whether
    or not the processor was an internal branch point.
    """
    def __init__(self, thisNewickStr, branchStr, options):
        Target.__init__(self)
        self.thisNewickStr = thisNewickStr
        self.options = options
        if self.options.seed != 'random':
            if branchStr == 'left':
                self.options.seed += 47
            elif branchStr == 'right':
                self.options.seed -= 61
            else:
                self.options.seed += 13
            self.options.seed = abs(self.options.seed)

    def run(self):
        #(name, step) = lsc.takeNewickStep( self.thisNewickStr, self.options )
        nt = newickTreeParser( self.thisNewickStr, 0.0 )
        name = lsc.nameTree(nt)
        commonParentDir = os.path.abspath(os.path.join( self.options.simDir, name ))
        lsc.myLog('TreeFollow.run(), %s %s\n' % (self.thisNewickStr, commonParentDir))
        if nt.distance == 0:            
            if nt.internal:
                lsc.myLog('TreeFollow.run() nt.distance == 0\n')
                # branch point
                branches = { 'left' : lsc.tree2str(nt.left),
                             'right': lsc.tree2str(nt.right) }
                for b in branches:
                    name = lsc.nameTree( newickTreeParser(branches[b], 0.0) )
                    childDir = os.path.join( self.options.simDir, name )
                    if not lsc.cycleIsLeaf( branches[b] ):
                        self.addChildTarget( Tree( branches[b], commonParentDir, b, self.options ))
        else:
            # stem with distance
            lsc.myLog('TreeFollow.run() nt.distance != 0.\n')
            self.addChildTarget( Tree( lsc.tree2str(nt), commonParentDir, 'stem', self.options ))
        self.addChildTarget( Transalign( commonParentDir, self.options ))
        self.addChildTarget( Stats( commonParentDir, self.options ))

class Cycle( Target ):
    """ The Cycle class creates the necessary directory structure for the
    given Cycle and then launches CycleStep1 as a child.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Target.__init__(self)
        self.thisNewickStr = thisNewickStr
        self.thisParentDir = thisParentDir
        self.thisStepSize  = thisStepSize
        self.thisDir = lsc.nameTree( newickTreeParser(thisNewickStr, 0.0) )
        self.options = options
    def run(self):
        os.mkdir( os.path.join( self.options.simDir, self.thisDir ))
        for d in ['inter', 'intra', 'logs', 'mobiles', 'stats']:
            os.mkdir( os.path.join( self.options.simDir, self.thisDir, d ))
        lsc.createNewCycleInfoXml( os.path.abspath(os.path.join(self.options.simDir, self.thisDir)), 
                                   self.thisParentDir, self.thisStepSize)
        self.addChildTarget( CycleStep1( self.thisNewickStr, self.thisParentDir, self.thisStepSize, self.options))

class CycleStep1( Cycle ):
    """ CycleStep1 consists of an evolver inter step and then the mobiles step.
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'CycleStep1' )
        cmd = [ lsc.which('touch') ]
        cmd.append( os.path.abspath(os.path.join( self.options.simDir, self.thisDir, 'cycleFollow1' )))
        lsc.runCommands( [cmd], 'CycleStep1' )
        
        
        
        self.setFollowOnTarget( CycleStep2( self.thisNewickStr, self.thisParentDir, 
                                            self.thisStepSize, self.options))

class CycleStep2( Cycle ):
    """ CycleStep2 
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'CycleStep2' )
        f = open( os.path.join( self.options.simDir, self.thisDir, 'cycleFollow2' ), 'w')
        f.write('words!\n')
        f.close()
        self.setFollowOnTarget( CycleStep3( self.thisNewickStr, self.thisParentDir, 
                                            self.thisStepSize, self.options))

class CycleStep3( Cycle ):
    """ CycleStep3 
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'CycleStep3' )
        f = open( os.path.join( self.options.simDir, self.thisDir, 'cycleFollow3' ), 'w')
        f.write('words!\n')
        f.close()
        self.setFollowOnTarget( CycleStep4( self.thisNewickStr, self.thisParentDir,
                                            self.thisStepSize, self.options))

class CycleStep4( Cycle ):
    """ CycleStep4 
    """
    def __init__(self, thisNewickStr, thisParentDir, thisStepSize, options):
        Cycle.__init__(self, thisNewickStr, thisParentDir, thisStepSize, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'CycleStep4' )
        f = open( os.path.join( self.options.simDir, self.thisDir, 'cycleFollow4' ), 'w')
        f.write('words!\n')
        f.close()

class Stats( Target ):
    """ The Stats object is a convenience class that launches
    StatsStep1 as a child.
    """
    def __init__(self, thisDir, options):
        Target.__init__(self)
        self.thisDir = thisDir
        self.options = options
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'Stats step' )
        f = open( os.path.join( self.thisDir, 'statsStep' ), 'w')
        f.write('statsStep\n')
        f.close()
        self.addChildTarget( StatsStep1( self.thisDir, self.options))

class StatsStep1( Stats ):
    """ StatsStep1 
    """
    def __init__(self, thisDir, options):
        Stats.__init__(self, thisDir, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'StatsStep1' )
        f = open( os.path.join( self.thisDir, 'statsStep1' ), 'w')
        f.write('statsStep\n')
        f.close()
        self.setFollowOnTarget( StatsStep2( self.thisDir, self.options))

class StatsStep2( Stats ):
    """ StatsStep2
    """
    def __init__(self, thisDir, options):
        Stats.__init__(self, thisDir, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'StatsStep2' )
        f = open( os.path.join( self.thisDir, 'statsStep2' ), 'w')
        f.write('statsStep\n')
        f.close()
        self.setFollowOnTarget( StatsStep3( self.thisDir, self.options))

class StatsStep3( Stats ):
    """ StatsStep3 
    """
    def __init__(self, thisDir, options):
        Stats.__init__(self, thisDir, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'StatsStep3' )
        f = open( os.path.join( self.thisDir, 'statsStep3' ), 'w')
        f.write('statsStep\n')
        f.close()
        self.setFollowOnTarget( StatsStep4( self.thisDir, self.options))

class StatsStep4( Stats ):
    """ StatsStep4 
    """
    def __init__(self, thisDir, options):
        Stats.__init__(self, thisDir, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'StatsStep4' )
        f = open( os.path.join( self.thisDir, 'statsStep4' ), 'w')
        f.write('statsStep\n')
        f.close()

class Transalign( Target ):
    """ The Transalign class is a convenience class that
    launches the two TransalignStep classes.
    """
    def __init__(self, thisDir, options):
        Target.__init__(self)
        self.thisDir = thisDir
        self.options = options
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'Transalign' )
        f = open( os.path.join( self.thisDir, 'transalignStep' ), 'w')
        f.write('transalignStep\n')
        f.close()
        self.addChildTarget( TransalignStep1( self.thisDir, self.options))

class TransalignStep1( Transalign ):
    """ TransalignStep1
    """
    def __init__(self, thisDir, options):
        Transalign.__init__(self, thisDir, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'TransalignStep1' )
        f = open( os.path.join( self.thisDir, 'transalignStep1' ), 'w')
        f.write('transalignStep1\n')
        f.close()
        self.setFollowOnTarget( TransalignStep2( self.thisDir, self.options))

class TransalignStep2( Transalign ):
    """ TransalignStep2
    """
    def __init__(self, thisDir, options):
        Transalign.__init__(self, thisDir, options)
    def run(self):
        lsc.verifyDirExists(os.path.abspath( os.path.join(self.options.simDir, self.thisDir )), 'TransalignStep2' )
        f = open( os.path.join( self.thisDir, 'transalignStep2' ), 'w')
        f.write('transalignStep2\n')
        f.close()
