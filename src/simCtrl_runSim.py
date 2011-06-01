#!/usr/bin/env python
"""
runSim.py
dent earl, dearl (a) soe ucsc edu
19 April 2009
In order to run a simulation you call runSim.py
runSim.py checks to make sure that every single
piece of software called (eventually) by the simulation
will exist and it does some prep work before the first
call to simTree.py, which begins the recursive process
of performing the genome simulation.
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
from evolverSimControl.lib.libSimControlClasses import SimTree
import evolverSimControl.lib.libSimControl as lsc
from jobTree.scriptTree.stack import Stack
from optparse import OptionParser
import os
from sonLib.bioio import newickTreeParser
import sys

# this is the first script to run in a simulation, so it will check for the
# existance of everything the entire simulation will end up calling, not
# just the scripts or external files used by runSim.py.
programs = ['cp',
            'mkdir',
            'evolver_evo', 'evolver_cvt', 'evolver_transalign',
            'evolver_drawrev', 'evolver_gff_cdsutr2exons.py',
            'evolver_gff_exons2introns.py', 'evolver_gff_featurestats2.sh',
            'evolver_gff_featurestats2.py', 'evolver_codon_report.pl',
            'evolver_merge_evostats.py', 'evolver_mobile_report.pl',
            'touch', 'ln', 'egrep', 'cat', 'simCtrl_completeTimestamp.py']
lsc.verifyPrograms(programs)
(CP_BIN, MKDIR_BIN) = programs[ 0:2 ]

def initOptions(parser):
    parser.add_option('--rootName',dest='rootName', 
                      help='name of the root genome, to differentiate it from the input newick.')
    # Sim Tree Options
    parser.add_option('-o', '--outDir',dest='outDir',
                      help='Out directory.')
    parser.add_option('-t', '--inputNewick',dest='inputNewick',
                      help='Newick tree.')
    parser.add_option('--testTree', action='store_true', 
                      default=False, dest='testTree',
                      help='Instead of performing a simulation, does dry run with empty dirs.')
    # Sim Control Options
    parser.add_option('--rootDir',dest='rootInputDir',
                      help='Input root directory.')
    parser.add_option('-s', '--step',dest='stepSize', action="store",
                      type ='float', default=0.001,
                      help='stepSize for each cycle. default=%default')
    parser.add_option('--params',dest='paramsDir',
                      help='Parameter directory.')
    parser.add_option('--seed',dest='seed',default='stochastic',
                      type='string', help='Random seed, either an int or "stochastic". default%default')
    parser.add_option('--noMEs', action='store_true', 
                      dest='noMEs', default=False, 
                      help=('Turns off all mobile element '
                            'and RPG modules in the sim. default=%default'))
    parser.add_option('--noGeneDeactivation', action='store_true', 
                      dest='noGeneDeactivation', default=False, 
                      help=('Turns off the gene deactivation step. '
                            'default=%default'))

def checkOptions(options, parser):
    if options.inputNewick is None:
        parser.error('Specify --inputNewick.')
    # check newickTree for reserved words
    nt = newickTreeParser(options.inputNewick, 0.0)
    if options.rootName is None:
        options.rootName = lsc.newickRootName(nt)
    else:
        options.rootName = options.rootName
    if newickContainsReservedWord(nt, options):
        parser.error('Newick tree contains reserved words.\n')
    
    # Sim Tree Options
    if options.outDir is None:
        parser.error('specify --outDir.\n')
    if os.path.exists(options.outDir):
       parser.error('%s already exists! If your simulation crashed, '
                    'relaunch it with "jobTreeRun --jobTree %s/" \n' % 
                    (os.path.join(options.outDir), options.jobTree))
    options.outDir = os.path.abspath(options.outDir)
    if not os.path.exists(options.outDir):
        os.mkdir(options.outDir)
    # Sim Control options
    if options.rootInputDir is None:
        parser.error('Specify --rootDir.\n')
    if not os.path.isdir(options.rootInputDir):
        parser.error('--rootDir "%s" not a directory!\n' % options.rootInputDir)
    options.rootInputDir = os.path.abspath(options.rootInputDir)
    
    if options.paramsDir is None:
        parser.error('Specify --params.\n')
    if not os.path.isdir(options.paramsDir):
        parser.error('Params dir "%s" not a directory!\n' % options.paramsDir)
    options.paramsDir = os.path.abspath(options.paramsDir)
    if options.stepSize <= 0:
        parser.error('specify positive stepSize.\n')
    if options.seed != 'stochastic':
        options.seed = int(options.seed)
        # otherwise we let the evolver tools choose their own seeds at random.

def newickContainsReservedWord(nt, options):
    """
    newickContainsReservedWord() checks the newick to make sure that
    there are no reserved names used as IDs. At present the only reserved
    name is 'parameters'.
    """
    reservedWords = set([ 'parameters', options.rootName ])
    if nt is None:
        return False
    if nt.iD in reservedWords:
        return True
    right = newickContainsReservedWord(nt.right, options)
    left  = newickContainsReservedWord(nt.left,  options)
    if left or right:
        return True
    return False

def checkForFiles(options):
    """ If files are missing, complains and dies.
    """
    if not os.path.exists(os.path.join(options.rootInputDir, 'seq.rev')):
        sys.stderr.write('Error, unable to find seq.rev in --rootDir %s.\n' % options.rootInputDir)
    for p in ['model.txt', 'model.mes.txt', 'mes.cfg']:
        if not os.path.exists(os.path.join(options.paramsDir, p)):
            sys.stderr.write('Error, unable to find %s in --params %s' % (p, options.paramsDir))
            sys.exit(1)

def populateRootDir(options):
    """ The first order of business in a simulation is to create the basic directory structure
    for the root genome and the parameters.
    """
    # mkdir is used here for simplicity in timing the creation of the diretory and 
    # subsequent two cp jobs for parameters.
    jobs = []
    jobs.append(['mkdir', '-p', os.path.join(options.outDir, 'parameters') ])
    lsc.runCommands(jobs, options.outDir)
    jobs = []
    jobs.append(['cp', '-r', options.rootInputDir, os.path.join(options.outDir, options.rootName)])
    jobs.append(['cp', os.path.join(options.paramsDir,'model.txt'), 
                 os.path.join(options.outDir, 'parameters')])
    jobs.append(['cp', os.path.join(options.paramsDir,'model.mes.txt'),
                 os.path.join(options.outDir, 'parameters')])
    if not options.noMEs:
        jobs.append(['cp', os.path.join(options.paramsDir,'mes.cfg'),
                     os.path.join(options.outDir, 'parameters')])
    lsc.runCommands(jobs, options.outDir, mode='p')
    options.paramsInputDir = options.paramsDir
    options.paramsDir    = os.path.abspath(os.path.join(options.outDir, 'parameters'))
    options.parentDir    = os.path.abspath(os.path.join(options.outDir, options.rootName))
    options.simDir, tail = os.path.split(options.parentDir)
    options.rootDir      = os.path.abspath(os.path.join(options.simDir, 'root'))
    lsc.createRootXmls(sys.argv, options)
    
def launchSimTree(options):
    jobResult = Stack(SimTree(options)).startJobTree(options)
    if jobResult:
        sys.stderr.write('Error, the jobTree contained %d failed jobs!\n' % jobResult)
        sys.exit(1)

def main():
    usage=('usage: %prog --rootName=name --parent=/path/to/dir --params=/path/to/dir '
           '--tree=newickTree --step=stepSize --out=/path/to/dir --seed=seedNumber\n\n'
           '%prog is used to initiate an evolver simulation using jobTree/scriptTree.')
    parser=OptionParser(usage=usage)
    initOptions(parser)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, parser)

    checkForFiles(options)
    populateRootDir(options)
    launchSimTree(options)

if __name__ == "__main__":
    main()
