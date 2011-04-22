#!/usr/bin/env python
"""
eval_evolverMAFcleaner.py
dent earl, dearl (a) soe ucsc edu
17 nov 2009
A script to clean the MAFs produced by evolver's CVT program.
Name.intra.chr31 becomes Name.c31
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
import os, re, sys

# def usage():
#     sys.stderr.write('USAGE: %s < file.maf > out.maf \n' % (sys.argv[0]))
#     sys.stderr.write('This script reads maf files from STDIN and writes\n'
#                      'to STDOUT. The names are changed from:\n'
#                      'Name.intra.chr31\n'
#                      'to:\n'
#                      'Name.c31\n')
#     sys.exit(2)

def main():

    # END OPTIONS
    ########################################
    #
    for line in sys.stdin:
        line=line.rstrip()
        line=line.replace('intra.chr', 'c')
        print('%s' %(line))

if __name__ == "__main__":
    main()
