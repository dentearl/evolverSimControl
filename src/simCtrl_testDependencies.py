#!/usr/bin/env python
"""
testDependencies.py
dent earl, dearl (a) soe ucsc edu
20 June 2011

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
import evolverSimControl.lib.libSimControl as lsc
import evolverSimControl.lib.libSimControlClasses as lscc
from jobTree.scriptTree.stack import Stack
from jobTree.scriptTree.target import Target
from optparse import OptionParser
import os
import re
from sonLib.bioio import newickTreeParser
import sys

def initOptions(parser):
    parser.add_option('--all', dest = 'all', default = False,
                      action = 'store_true',
                      help = 'Runs all checks. default=%default')
    parser.add_option('--mafJoin', dest = 'mafJoin', default = False,
                      action = 'store_true',
                      help = 'Checks for mafJoin. default=%default')

def testSimple():
    lsc.verifyPrograms(lsc.requiredPrograms)

def testMafJoin():
    programs = ['mafJoin']
    lsc.verifyPrograms(programs)

def checkOptions(options, parser):
    pass

def main():
    usage=('usage: %prog [options]\n\n'
           '%prog checks to ensure that all required executables are contained in the PATH.')
    parser = OptionParser(usage = usage)
    initOptions(parser)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, parser)

    testSimple()
    print 'Dependencies for simple simulations..... OK'
    
    if options.mafJoin or options.all:
        testMafJoin()
        print 'Dependencies for maf construction....... OK'

if __name__ == '__main__':
    main()
