#!/usr/bin/env python
"""
simCtrl_postSimReAligner.py
dent earl, dearl (a) soe ucsc edu
12 april 2010

This script will correct for a mistake in simulation
setup. IF your simulation began with a stem, then it will be
impossible to infer the root genome; the best inference would
lead you to the branch point before the root, but not up the
stem. This script allows you to automatically generate
replacement alignments from leaves to root-branch via the
leaves to root alignment.
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
from optparse import OptionParser
import os
import sys
import subprocess
import simulation.lib.libSimControl as LSC

programs = ['evolver_transalign']
LSC.verifyPrograms(programs)
TRANS_BIN = programs[0]

def usage():
    print 'USAGE: '+sys.argv[0]+\
          ' --newRoot [root-branch.aln] --logDir [dir] '+\
          ' leaf1.aln leaf2.aln leaf3.aln ...'
    print __doc__
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-l', '--logDir',dest='logDir', default=os.getcwd(),
                      help='Output directory, only for logfiles.')
    parser.add_option('-n', '--newRoot',dest='newRoot',
                      help='The new root .aln.')
    parser.add_option('-d', '--debug', dest='debug', default=False,
                      action='store_true', help='Prints command, but does not execute')

def checkOptions(options, args):
    if not options.newRoot:
        sys.stderr.write('%s: Error, you must specify --newRoot' %(sys.argv[0]))
        usage()
    if options.logDir:
        if not os.path.isdir(options.logDir):
            sys.stderr.write('%s: Error, directory "%s" is not a directory!\n' % (sys.argv[0], options.logDir))
            usage()
    else:
        options.logDir = os.path.abspath(os.getcwd())
    options.newRoot = os.path.abspath(options.newRoot)
    options.logDir  = os.path.abspath(options.logDir)
    for a in args:
        if not os.path.exists(a):
            sys.stderr.write('%s: Error, alignment "%s" does not exist!\n' % (sys.argv[0], a))
            usage()

def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options, args)
    for a in args:
        homeDir=os.path.dirname(a)
        genome=os.path.basename(homeDir)
        CMD=[]
        CMD.append(TRANS_BIN)
        CMD.append('-in1')
        CMD.append('%s' %(options.newRoot))
        CMD.append('-in2')
        CMD.append('%s' %(a))
        CMD.append('-out')
        CMD.append('%s' %(os.path.join(homeDir,'newRoot.aln.rev')))
        CMD.append('-log')
        CMD.append('%s' %(os.path.join(options.logDir, 'newRoot.'+genome+'.aln.log')))
        if options.debug:
            print CMD
        else:
            subprocess.Popen(CMD)

if __name__ == "__main__":
    main()
