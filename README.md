#evolverSimControl
(c) 2009 - 2011 The Authors, see LICENSE.txt for details.

##Authors
[Dent Earl](https://github.com/dentearl/), [Benedict Paten](https://github.com/dentearl/), Mark Diekhans

The evolver team is responsible for items in external/ : George Asimenos and [Robert C. Edgar](http://www.drive5.com/), Serafim Batzoglou and [Arend Sidow](http://mendel.stanford.edu/sidowlab/).

## Summary
A [jobTree](https://github.com/benedictpaten/jobTree/) based simulation manager for the [Evolver](http://www.drive5.com/evolver/) genome evolution simulation tool suite. 

**evolverSimControl** (eSC) can be used to simulate multi-chromosome genome evolution on an arbitrary phylogeny ([Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html)). In addition to simply running evolver, **eSC** also automatically creates statistical summaries of the simulation as it runs including text and image files. Also included are convenience scripts to: check on a running simulation and see detailed status and logging information; extract fasta sequence files from the leaf nodes of a completed simulation; extract pairwise multiple alignment files ([.maf](http://genome.ucsc.edu/FAQ/FAQformat.html#format5)) from leaf and branch nodes from a completed simulation and with the help of [mafJoin](https://github.com/dentearl/mafTools/), join them together into a single maf covering the entire simulation.

The use of jobTree means that you can run **eSC** on a cluster running a jobTree supported batch system, on a multi-cored server or on your laptop.

##Dependencies
* **sonLib**: https://github.com/benedictpaten/sonLib/
* **jobTree**: https://github.com/benedictpaten/jobTree/
* **evolver**: http://www.drive5.com/evolver/ Specifically, make sure that the Evolver tools are on your <code>PATH</code> environmental variable and that their names are preceeded with <code>evolver_</code>. Specifically all of the following list of files need to be on your <code>PATH</code>.
    * <code>evolver_cvt</code>
    * <code>evolver_evo</code>
    * <code>evolver_transalign</code>
* **trf**: http://tandem.bu.edu/trf/trf.html Tandem Repeats Finder.
* **mafJoin**: https://github.com/dentearl/mafTools Not necessary for simple simulations, mafJoin (part of mafTools) is only needed if you wish to create a maf alignment of all sequences following a simulation.
* **R**: http://cran.r-project.org/ Only necessary if you wish to use the <code>simCtrl_postSimAnnotDistExtractor.py</code> script to view annotation size distributions following a simulation.
* **ggplot2** for R: in R type <code>install.packages("ggplot2")</code> Only necessary if you wish to use the <code>simCtrl_postSimAnnotDistExtractor.py</code> script to view annotation size distributions following a simulation.

##Requirements
* Linux on i86 Intel. This is due to core Evolver executables being distributed as pre-compiled binaries.

##Installation
1. Download the package. Consider making it a sibling directory to <code>jobTree/</code> and <code>sonLib/</code>.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.
4. Edit your <code>PYTHONPATH</code> environmental variable to contain the parent directory of the <code>evolverSimControl/</code> directory.
5. Type <code>make test</code>.

##Example
This example will work you through a small simulation using the toy test example available at http://soe.ucsc.edu/~dearl/software/evolverSimControl/. If you want to create your own infile you can use [evolverInfileGeneration](https://github.com/dentearl/evolverInfileGeneration) to generate your own infile set.

1. Download and expand the toy archive. For simplicity I'll assume that both <code>root/</code> and <code>params/</code> are in the working directory, i.e. <code>./</code> .
2. Next we run the runSim program:
    * <code>$ simCtrl_runSim.py --inputNewick '(Knife:0.004, (Fork:0.003, (Ladle:0.002, (Spoon:0.001, Teaspoon:0.001)S-TS:.001)S-TS-L:.001)S-TS-L-F:0.001);' --outDir toyExampleSim --rootDir root/ --rootName hg18 --paramsDir params/ --jobTree jobTreeToyExampleSim --maxThreads 32 --seed 3571</code>
    * You can check on a running simulation by using <code>simCtrl_checkSimStatus.py</code> , use <code>--help</code> for options.
3. Post simulation you can run <code>simCtrl_postSimFastaExtractor.py</code> to extract fasta sequence files from the genomes.
4. You may also wish to run <code>simCtrl_postSimAnnotDistExtractor.py</code> which will use the ggplot2 package for R to display the length distributions of some of the annotations.
5. You may also wish to construct a single maf for the simulation using <code>simCtrl_postSimMafExtractor.py</code> which will use [mafJoin](https://github.com/dentearl/mafTools/) to join the pairwise maf output from Evolver into a single simulation wide maf. This process is extremely memory intensive with the 120Mb Mammal simulation eventually requiring aprroximately 250Gb of memory.

##Use
###Initiating a simulation
In order to run eSC you will need an infile set, a parameter set, a phylogenetic tree and optionally a mobile element library and mobile element parameter set. Infile sets can be created using [evolverInfileGenerator](https://github.com/dentearl/evolverInfileGenerator/) or from scratch. Parameter sets can be generated by reading primary literature and coming up with reasonable values. Phylogenetic trees need to be in Newick format.

Available options for running a simulation are listed below.

<code>$ bin/simCtrl_runSim.py --help</code>

<code>Usage: simCtrl_runSim.py --rootName=name --rootDir=/path/to/dir --paramsDir=/path/to/dir --tree=newickTree --stepLength=stepLength --outDir=/path/to/dir --jobTree=/path/to/dir [options]</code>

simCtrl_runSim.py is used to initiate an evolver simulation using jobTree/scriptTree.

Options:

* <code>-h, --help</code> show this help message and exit
* <code>--rootDir=ROOTINPUTDIR</code> Input root directory.
* <code>--rootName=ROOTNAME</code> name of the root genome, to differentiate it from the input Newick. default=root
* <code>--inputNewick=INPUTNEWICK</code> Newick tree. http://evolution.genetics.washington.edu/phylip/newicktree.html
* <code>--stepLength=STEPLENGTH</code> stepLength for each cycle. default=0.001
* <code>--paramsDir=PARAMSDIR</code> Parameter directory.
* <code>--outDir=OUTDIR</code> Out directory.
* <code>--seed=SEED</code> Random seed, either an int or "stochastic". default=stochastic
* <code>--noMEs</code> Turns off all mobile element and RPG modules in the sim. default=False
* <code>--noBurninMerge</code> Turns off checks for an aln.rev file in the root dir. default=False
* <code>--noGeneDeactivation</code> Turns off the gene deactivation step. default=False
* <code>--maxThreads=MAXTHREADS</code> The maximum number of threads to use when running in single machine mode. default=4
* ... and all other jobTree standard options.

###Simulation Status
To check on a running simulation you can use the <code>simCtrl_checkSimStatus.py</code> script.

<code>$ bin/simCtrl_checkSimStatus.py --help</code>

<code>Usage: simCtrl_checkSimStatus.py --simDir path/to/dir [options]</code>

simCtrl_checkSimStatus.py can be used to check on the status of a running or completed
evolverSimControl simulation.

Options:

* <code>-h, --help</code> show this help message and exit
* <code>--simDir=SIMDIR</code> Parent directory.
* <code>--drawText, --drawTree</code> prints an ASCII representation of the current tree status. default=False
* <code>--curCycles</code> prints out the list of currently running cycles. default=False
* <code>--stats</code> prints out the statistics for cycle steps. default=False
* <code>--cycleStem</code> prints out a stem and leaf plot for completed cycle runtimes, in seconds. default=False
* <code>--cycleStemHours</code> prints out a stem and leaf plot for completed cycle runtimes, in hours. default=False
* <code>--printChrTimes</code> prints a table of chromosome lengths (bp) and times (sec) for intra chromosome evolution step (CycleStep2).
* <code>--cycleList</code> prints out a list of all completed cycle runtimes. default=False
* <code>--html</code> prints output in HTML format for use as a cgi. default=False
* <code>--htmlDir=HTMLDIR</code> prefix for html links.

###Sequence Extraction
To extract fasta sequences from a completed simulation you can use the <code>simCtrl_postSimFastaExtractor.py</code> script.

<code>$ bin/simCtrl_postSimFastaExtractor.py --help</code>

<code>Usage: simCtrl_postSimFastaExtractor.py --simDir path/to/dir [options]</code>

simCtrl_postSimFastaExtractor.py takes in a simulation directory and then extracts the sequences
of leaf nodes in fasta format and stores them in the respective step's directory.

Options:

* <code>-h, --help</code> show this help message and exit
* <code>--simDir=SIMDIR</code> the simulation directory.
* <code>--allCycles</code> extract fastas from all cycles, not just leafs. default=False

###Simulation maf creation
To create a single maf reflecting the evolutionary history of the entire simulation <code>simCtrl_postSimFastaExtractor.py</code> script.

<code>$ bin/simCtrl_postSimMafExtractor.py --help</code>

<code>Usage: simCtrl_postSimMafExtractor.py --simDir path/to/dir [options]</code>

simCtrl_postSimMafExtractor.py requires mafJoin which is part of mafTools and is available
at https://github.com/dentearl/mafTools/ . 

Options:

* <code>-h, --help</code> show this help message and exit
* <code>--simDir=SIMDIR</code> Simulation directory.
* <code>--maxBlkWidth=MAXBLKWIDTH</code> Maximum mafJoin maf block output size. May be reduced towards 250 for complicated phylogenies. default=10000
* <code>--maxInputBlkWidth=MAXINPUTBLKWIDTH</code> Maximum mafJoin maf block input size. mafJoin will cut inputs to size, may result in long
runs for very simple joins. May be reduced towards 250 for complicated phylogenies.
default=1000                        
* <code>--noBurninMerge</code> Will not perform a final merge of simulation to the burnin. default=False
* <code>--maxThreads=MAXTHREADS</code> The maximum number of threads to use when running in single machine mode. default=4
* ... and all other jobTree standard options.
