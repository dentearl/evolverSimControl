#!/usr/bin/env python
"""
eval_evolverMAFcleaner.py
dent earl, dearl (a) soe ucsc edu
17 nov 2009
A script to clean the MAFs produced by evolver's CVT program.
Name.intra.chr31 becomes Name.c31
"""
##############################
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
