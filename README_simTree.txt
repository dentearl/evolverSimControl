README simTree.txt
dent earl, dearl (a) soe ucsc edu
20 Oct 2009

Documentation for the simTreeXXX.py package of scripts.

Some of these scripts must be run in the jobTree.py framework.
For a bash-python set of scripts that does all of this 
in on a single machine, instead of a cluster, see
simTreeStandAlone.py and cycleStandAlone.sh.

The simTreeXXX.py set of scripts constitute a recursive
suite for walking down a newick tree and running one
evolver cycle at each step, two cycles if at a branch
point.

##############################
Members
##############################
 simTree.py
 simTreeFollowUp.py
 simTreeStandAlone.py

##############################
Descriptions
##############################
# simTree.py
##############################
simTree.py is intended to be an end-user interface for
the cycleXXX.py scripts (and therefore the 'evolver' 
genome simulation suite of tools.)
simTree.py is the first of two scripts that automatically
simulations an entire phylogenetic tree. Simulations must
start with a root directory (which must contain both an
annotation file (called annots.gff) and a sequence file
(called seq.rev). Also required is a parameters directory
which needs to contain a file called mes.fa (moble element
library in fasta format) and model.txt (evolver parameter
file).
simTree.py also requires a newick tree from which to simulate.
simTree.py creates up to two jobTree child processes, either
one call to cycleBegin.py if descending down a single branch,
or two calls to cycleBegin.py if at a branch point in the tree.
There is also an --out option to specify the directory in which
all child cylces will be placed. The default is to place the
child cycles in the same directory as the root genome.
simTree.py has one jobTree follow up job, which is a call to
simTreeFollowUp.py. It is important to note that simTree.py
reads the newick tree and passes it on, unchanged, to
simTreeFollowUp.py
Example run (reference only):
./simTree.py  --parent parentDir/ \
              --params paramsDir/ \
              --tree '(Steve:0.003, (Zack:.005, Chris:.011):.001):.002;' \
              --step 0.02 \
              --out resultsDirectory/ \
              --jobFile JOB_FILE
##############################
# simTreeFollowUp.py
##############################
simTreeFollowUp.py takes the same inputs as simTree.py and alters
the inputed newick tree, cutting off whatever step was taken by 
simTree.py. 
simTreeFollowUp.py creates up to two jobTree child processes, either
one call to simTree.py if descending down a single branch,
or two calls to simTree.py if at a branch point in the tree.
simTreeFollowUp.py has no follow up jobTree commands.

Example run (reference only):
./simTreeFollowUp.py --parent parentDir/ \
                     --params paramsDir/ \
                     --tree '(Steve:0.003, (Zack:.005, Chris:.011):.001):.002;' \
                     --step 0.02 \
                     --out resultsDirectory/ \
                     --jobFile JOB_FILE
##############################
# simTreeStandAlone.py
##############################
simTreeStandAlone.py replicates the functionality of simTree.py, but does so
without using the jobTree.py framework. simTreeStandAlone.py is suitable for
running on lab size servers (8+ cores). It is a recursive script that issues 
backgrounded subprocesses to a Bash script that actually runs the cycle 
(cycleStandAlone.sh). All branches are run in parallel.
It takes all the same input as simTree.py except there is no --jobFile input.
Additionally, simTreeStandAlone.py will write to a file in the 'parent'
directory (provided it is not the root) the next recursive command(s) that
are to be run. These commands will either be in files named 'nextCommand.txt'
in the case of descending down a steam, or 'nextCommand_left.txt' and
'nextCommand_right.txt' in the case of a branch point. Should there be a 
downstream failure, the user should be able to take these commands and 
reissue the job using the 'parent' directory as the new root node.
Example run (reference only):
./simTreeStandAlone.py  --parent parentDir/ \
                        --params paramsDir/ \
                        --tree '(Steve:0.003, (Zack:.005, Chris:.011):.001):.002;' \
                        --step 0.02 \
                        --out resultsDirectory/ \
