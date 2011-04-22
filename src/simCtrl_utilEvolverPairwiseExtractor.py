#!/usr/bin/env python
"""
eval_evolverPairwiseExtractor.py
dent earl, dearl (a) soe ucsc edu
16 nov 2009
A script that will take the path to a simCtrl_simTree.py
out directory and make repeated calls to transalign
and to cvt (in fabulous PARALLEL-o-vision!) in order
to calculate all leaf pairwise alignments and generate
all pairwise MAFs (or Fastas...).
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
import glob, os, re, subprocess, sys
import simulation.lib.libSimControl as LSC

programs = ['evolver_cvt', 'evolver_transalign', 'eval_evolverMAFcleaner.py']
LSC.verifyPrograms(programs)
(CVT_BIN, TRANS_BIN, MAFCLEANER) = programs

def usage():
    sys.stderr.write('USAGE: %s --cycleDir <dir> --out <dir> [optional: --leavesOnly --log]\n' % (sys.argv[0]))
    print __doc__
    sys.exit(2)

def logMessage(message):
    curr = os.curdir
    logPath = os.path.join(curr, 'pairwise_log.log')
    if(os.path.exists(logPath)):
        os.remove(logPath)
    FILE = open(logPath, 'w')
    FILE.write( '%s' % (message))
    FILE.close()
        
def main():
    parser=OptionParser()
    parser.add_option('-i', '--cycleDir',dest='inDir',
                      help='Out directory from simTree.py.')
    parser.add_option('-o', '--out',dest='outDir',
                      help='Location to place all pairwise aln.revs, MAFs.')
    parser.add_option('-l', '--log', action='store_true', dest='log',
                      default=False, help='Record output from evolver suite information.')
    parser.add_option('-s', '--leavesOnly', action='store_true', dest='leavesOnly',
                      default=False, help='Will exclude any directory with a "-" in the name, which are taken to be internal nodes.')
    parser.add_option('-m', '--maf', action='store_true', dest='isMAF',
                      default=True, help='Extracts pairwise MAF.')
    (options, args) = parser.parse_args()

    if (options.inDir == None):
        sys.stderr.write('%s: Error, specify --cycleDir.\n' % sys.argv[0])
        usage()
    if (options.outDir == None):
        options.outDir = options.inDir
    if (not os.path.isdir(options.inDir)):
        sys.stderr.write('%s: Error, cycle directory "%s" not a directory!\n' % (sys.argv[0], options.inDir))
        usage()
    if(options.outDir != None):
        if(not os.path.exists(options.outDir)):
            os.mkdir(options.outDir)
    outDir = os.path.abspath(options.outDir)
    inDir = os.path.abspath(options.inDir)
    # END OPTIONS
    ########################################
    # 
    cycles = glob.glob(inDir+'/*')
    cmdArray=[]
    jobsArray=[]
    knownPairs={}
    dash = re.compile('-')
    #####
    # MAIN Extraction
    for i in cycles:
        if(not os.path.isdir(i)):
                continue
        for j in cycles:
            if(not os.path.isdir(j)):
                continue
            if (i != j):
                if (i+' - '+j not in knownPairs) and (j+' - '+i not in knownPairs):
                    if(options.leavesOnly):
                        if(dash.search(i) or dash.search(j)):
                            continue
                    knownPairs[i+' - '+j]=1
                    nameA    = str(i.split('/')[-1:])
                    nameB    = str(j.split('/')[-1:])
                    nameA    = nameA.replace('[','')
                    nameA    = nameA.replace(']','')
                    nameA    = nameA.replace('\'','')
                    nameB    = nameB.replace('[','')
                    nameB    = nameB.replace(']','')
                    nameB    = nameB.replace('\'','')
                    outName  = nameA+'.'+nameB+'.aln.rev'
                    outTrans = os.path.join(outDir,outName)
                    outName  = nameA+'.'+nameB+'.aln.maf'
                    outNameC = nameA+'.'+nameB+'.aln.clean.maf'
                        
                    outCVT   = os.path.join(outDir,outName)
                    outCLEAN = os.path.join(outDir,outNameC)
                    cmd = TRANSALIGN +\
                          ' -in1 '+i+'/root.aln.rev'+\
                          ' -in2 '+j+'/root.aln.rev'+\
                          ' -out '+outTrans+'; '+\
                          CVT +\
                          ' -fromrev '+outTrans+\
                          ' -tomaf '+outCVT+\
                          '; '+\
                          MAFCLEANER+\
                          ' < '+outCVT+\
                          ' > '+outCLEAN
                    
                    cmdArray.append(cmd)
                    if (options.log):
                        curr = os.curdir
                        logPath = os.path.join(curr, 'pairwise_log.log')
                    else:
                        logPath = os.path.abspath('/dev/null')
                    logFILE = open(logPath, 'w')
                    job = subprocess.Popen(cmd, shell=True, stderr=logFILE, stdout=logFILE)
                    jobsArray.append(job)
    i=-1
    for job in jobsArray:
        job.wait()
        i=i+1
        if(job.returncode):
            if(options.log):
                logMessage('%s: Error in a sub-processes "%s", returncode: %s.\n' %(sys.argv[0], cmdArray[i], job.returncode))
            else:
                sys.stderr.write('%s: Error in a sub-processes "%s", returncode: %s.\n' %(sys.argv[0], cmdArray[i], job.returncode))
            sys.exit(1)

if __name__ == "__main__":
    main()
