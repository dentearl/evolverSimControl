#!/usr/bin/env python
# eval_maskedFaLength.py
# dent earl, dearl (a) soe ucsc edu
# 16 dec 2009
#
# scan through a FASTA, counting up the 
# unmasked sequence (i.e. CAPITAL ACTG)
#
##############################
import os, re, sys

def usage():
    sys.stderr.write('USAGE: %s file.fa\n' %(sys.argv[0]))
    sys.exit(2)

def main():
    if(len(sys.argv) != 2):
        usage()
    if(not os.path.isfile(sys.argv[1])):
        usage()
    infile = open(sys.argv[1], 'r')
    patCmnt= re.compile('^>')
    DNA = re.compile('[ACTG]')
    fastaLen = 0
    while infile:
        line = infile.readline()
        if(not line):
            break
        if(patCmnt.match(line)):
            continue
        fastaLen = fastaLen + len(DNA.findall(line))

    print fastaLen

if __name__ == "__main__":
    main()
