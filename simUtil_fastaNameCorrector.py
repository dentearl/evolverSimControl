#!/usr/bin/env python
"""
eval_fastaNameCorrector.py
dent earl, dearl (a) soe ucsc edu
17 nov 2009
A script that reads fasta files from STDIN and writes
to STDOUT. The comment line is changed from:
>chr31
to:
>name.c31
"""
##############################
from optparse import OptionParser
import os, re, sys

def usage():
    sys.stderr.write('USAGE: %s --name <name> < file.fa > out.fa \n' % (sys.argv[0]))
    sys.stderr.write('This script reads fasta files from STDIN and writes\n'
                     'to STDOUT. The comment line is changed from:\n'
                     '>chr31\n'
                     'to:\n'
                     '>[nameOption].c31\n')
    sys.exit(2)

def main():
    parser=OptionParser()
    parser.add_option('-i', '--name',dest='name',
                      help='The name to be appended at the start of the comment line.')
    (options, args) = parser.parse_args()
    if (options.name == None):
        sys.stderr.write('%s: Error, specify --name.\n' % sys.argv[0])
        usage()
    # END OPTIONS
    ########################################
    #
    for line in sys.stdin:
        line=line.rstrip()
        line=line.replace('>', '>'+options.name+'.')
        print('%s' %(line))

if __name__ == "__main__":
    main()
