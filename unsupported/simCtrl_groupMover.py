#!/usr/bin/env python
"""
groupMover.py
 dent earl, dearl (a) soe ucsc edu

groupMove.py is a python wrapper for the evolver suite of genome
evolution tools. groupMove takes a destination directory and then
moves every file passed as an argument to that directory. It is
intended to be used to move files from a local temp directory on a
cluster node to the hive filesystem after some cluster job.

2 March 2010
"""
######################################
import os, sys
from optparse import OptionParser

def usage():
    print "USAGE: %s --dest [destination directory] file1 file2 file3 ..." %(sys.argv[0])
    print __doc__
    sys.exit(2)

def main(argv):
    parser=OptionParser()
    parser.add_option('-d', '--dest',dest='destDir',
                      help='Destination directory.')
    (options, args) = parser.parse_args()
    if (options.destDir == None):
        sys.stderr.write('%s: Error, specify destination dir.\n' % sys.argv[0])
        usage()
    if (not os.path.isdir(options.destDir)):
        sys.stderr.write('%s: Error, Destination dir "%s" not a directory!\n' % (sys.argv[0], options.destDir))
        usage()
    for a in args:
        a = os.path.abspath(a)
        if(not os.path.exists(a)):
            sys.stderr.write('%s: Error, file to move, %s, does not exist.\n' % (sys.argv[0], a))
        os.rename(a, os.path.join(options.destDir, os.path.basename(a)))

if __name__ == '__main__':
    main(sys.argv)
