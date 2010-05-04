#!/usr/bin/env python
"""
simCtrl_wrapperTRF.py, 17 February 2009
 dent earl, dearl (a) soe ucsc edu

simCtrl_wrapperTRF.py takes all the the same arguments that
trf takes and passes them to trf. wrapper_TRF then checks the
return code and acts appropriately.

trf has issues with return codes and this is a simple way around
some of the issues.

the return code for trf  is the number of seqs it masked, or 255
if there was an error, [-1 is turned into an unsigned int] (if
you tried to mask 255 seqs you're out of luck. We only do one
at a time, so we expect a retcode of 1.)

"""
########################################
import glob, os, subprocess, sys
import simulation.lib.libSimControl as LSC

programs = ['trf']
LSC.verifyPrograms(programs)
TRF_BIN = programs[0]

def usage():
    print "USAGE: %s [trf options]" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def main(argv):
    cmd=[TRF_BIN]
    cmd.extend(argv[1:])
    localDir=os.path.dirname(argv[1])
#    sys.stderr.write('I want to run:%s\nfrom this directory:%s\n' %(cmd, localDir))
#    myDir=glob.glob(os.path.join(localDir,'*'))
#    sys.stderr.write('contents of that directory: %s\n' %(myDir))
#    subprocess.Popen('ls -al '+localDir, shell=True)
    p = subprocess.Popen(cmd, cwd=localDir)
    p.wait()
    if p.returncode != 1:
        if p.returncode < 0:
            sys.stderr.write('%s: Experienced an error while trying to execute: %s SIGNAL:%d\n' %(sys.argv[0], ' '.join(cmd), -(p.returncode)))
        else:
            sys.stderr.write('%s: Experienced an error while trying to execute: %s retcode:%d\n' %(sys.argv[0], ' '.join(cmd), p.returncode))
        sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
    main(sys.argv)
