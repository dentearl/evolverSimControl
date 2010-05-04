cycle.txt
dent earl, dearl (a) soe ucsc edu
20 Oct 2009

Documentation for the cycleXXXX.py package of scripts.

These scripts must be run in the jobTree.py framework.
For a bash-python set of scripts that does all of this
in on a single machine, instead of a cluster, see
cycleStandAlone.sh and simTreeStandAlone.py.

These scripts are automatically called by simTree.py,
this documentation is provided as a reference.

##############################
Members
##############################
 cycleBegin.py
 cycleMid_1.py
 cycleMid_2.py
 cycleEnd.py

##############################
Descriptions
##############################
 cycleBegin.py
##########
cycleBegin.py starts a cycle. It requires a parent directory that
contains a file called annots.gff and a file called seq.rev . It 
also requires a parameters directory which contains a file called
mes.fa (moble element library in fasta format) and model.txt 
(evolver parameter file). You must specify a child which will become
a new cycle directory. You may specify the step size of the cycle,
default is 0.001. You must specify the jobFile (jobTree.py formatted
xml file).
cycleBegin will call evo64 for the initial interchromosomal evolution
simulation step. The follow up job is cycleMid_1.py
Example call:
cycleBegin.py --parent rootDir/ --child childDir/\
              --params parameterDir --step 0.002 --jobFile JOB_FILE

 cycleMid_1.py
##########
cycleMid_1.py takes all of the same inputs as cycleBegin.py.
It creates the following children:
 1) transalign the root-to-parent alignment with the parent-to-inter
    alignment in order to get a root-to-inter alignment.
 2) for every chromosome, start an intrachromosomal evolution sim.
cycleMid_1.py has a follow up job, cycleMid_2.py
Example call:
cycleMid_1.py --parent rootDir/ --child childDir/\
              --params parameterDir --step 0.002 --jobFile JOB_FILE

 cycleMid_2.py
##########
cycleMid_2.py takes all of the same inputs as cycleMid_1.py, except
for --step, which is not accepted.
It creates the following children:
 1) concatenate all chromosome .gff files into one .gff file.
 2) merge all chromosome alignment files into one alignment file.
 3) merge all chromosome sequence files into one sequence file.
cycleMid_2.py has a follow up job, cycleEnd.py
Example call:
cycleMid_2.py --parent rootDir/ --child childDir/\
              --params parameterDir --jobFile JOB_FILE

 cycleEnd.py
##########
cycleEnd.py takes all of the same inputs as cycleMid_2.py.
It creates the following children:
 1) transalign root-to-inter alignment with inter-to-intra alignment
    in order to get a root-to-child alignment (input for descendant
    cycles).
 2) align cds sequences.
cycleEnd.py has no follow up jobs.
cycleEnd.py --parent rootDir/ --child childDir/\
              --params parameterDir --jobFile JOB_FILE
