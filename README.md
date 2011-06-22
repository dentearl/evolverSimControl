#evolverSimControl
(c) 2009 - 2011 The Authors, see LICENSE.txt for details.

##Authors
[Dent Earl](https://github.com/dentearl/), [Benedict Paten](https://github.com/dentearl/), Mark Diekhans

The evolver team is responsible for items in external/ : George Asimenos, Robert C. Edgar, Serafim Batzoglou and Arend Sidow.

## Summary
A [jobTree](https://github.com/benedictpaten/jobTree/) based wrapper for the [Evolver](http://www.drive5.com/evolver/) genome evolution simulation tool suite.

##Dependencies
* sonLib: https://github.com/benedictpaten/sonLib/
* jobTree: https://github.com/benedictpaten/jobTree/
* Evolver: http://www.drive5.com/evolver/ Specifically, make sure that the Evolver tools are on your <code>PATH</code> environmental variable and that their names are preceeded with <code>evolver_</code>. Specifically all of the following list of files need to be on your <code>PATH</code>. Depending on how you add these things to your <code>PATH</code>, you may need to edit some of the evolver python scripts to ensure that they can import the evolver <code>gff.py</code> library.
    * <code>evolver_codon_report.pl</code> *
    * <code>evolver_cvt</code>
    * <code>evolver_drawrev</code> *
    * <code>evolver_evo</code>
    * <code>evolver_evostats_report.py</code>
    * <code>evolver_gff_cdsutr2exons.py</code>
    * <code>evolver_gff_exons2introns.py</code>
    * <code>evolver_gff_featurestats2.py</code>
    * <code>evolver_gff_featurestats2.sh</code> *
    * <code>evolver_gene_deactivate.sh</code> *
    * <code>evolver_handle_mobiles.pl</code> *
    * <code>evolver_merge_evostats.py</code>
    * <code>evolver_mobile_report.pl</code> *
    * <code>evolver_transalign</code>
    * <code>evolver_trf2gff.py</code>
    * _* indicates this script/bin is included in the <code>external/</code> directory._
* trf: http://tandem.bu.edu/trf/trf.html Tandem Repeats Finder.
* mafJoin: https://github.com/dentearl/mafTools Not necessary for simple simulations, mafJoin (part of mafTools) is only needed if you wish to create an maf alignment of all sequences following a simulation.
* R: http://cran.r-project.org/ Only necessary if you wish to use the <code>simCtrl_postSimAnnotDistExtractor.py</code> script to view annotation size distributions following a simulation.
* ggplot2 for R: in R type <code>install.packages("ggplot2')</code> Only necessary if you wish to use the <code>simCtrl_postSimAnnotDistExtractor.py</code> script to view annotation size distributions following a simulation.

##Requirements
* Linux on i86 Intel. This is due to core Evolver executables being distributed as precompiled binaries.

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
    * <code>$ simCtrl_runSim.py --inputNewick '(Knife:0.004, (Fork:0.003, (Ladle:0.002, (Spoon:0.001, Teaspoon:0.001)S-TS:.001)S-TS-L:.001)S-TS-L-F:0.001);' --outDir toyExampleSim --rootDir root/ --rootName hg18 --params params/ --jobTree jobTreeToyExampleSim --maxThreads 32 --seed 3571</code>
    * <code>--outDir</code> is where you simulation is going to end up.
    * <code>--rootDir</code> should point to the <code>root/</code> dir you created.
    * <code>--rootName</code> in this case is hg18. It's set in the infile creation step and you can pull this out of a .rev file with <code>evolver_cvt -dumpchrids path/to/seq.rev</code>
    * <code>--params</code> should point to the <code>params/</code> dir you created.
    * <code>--noMEs</code> turns off mobile element library simulation. If you leave this **out** then the <code>params/</code> dir must contain <code>mes.cfg</code> and <code>model.mes.txt</code>, and the <code>root/</code> dir must contain a directory named <code>mobiles/</code> that contains the files <code>LTR.fa</code>, <code>ME.fa</code>, and <code>ME.gff</code>.
    * <code>--maxThreads</code> is a [jobTree](https://github.com/benedictpaten/jobTree/) option for limiting the maximum number of parallel threads. The default is rather low.
    * <code>--seed</code> allows you to give a random seed (an integer) to the simulation. The default is the string 'stochastic'.
    * You can check on a running simulation by using <code>simCtrl_checkSimStatus.py</code> , use <code>--help</code> for options.
3. Post simulation you can run <code>simCtrl_postSimFastaExtractor.py</code> to extract fasta sequence files from the genomes.
4. You may also wish to run <code>simCtrl_postSimAnnotDistExtractor.py</code> which will use the ggplot2 package for R to display the length distributions of some of the annotations.
5. You may also wish to construct a single multiple alignment file ([.maf](http://genome.ucsc.edu/FAQ/FAQformat.html#format5)) for the simulation using <code>simCtrl_postSimMafExtractor.py</code> which will use [mafJoin](https://github.com/dentearl/mafTools/) to join the pairwise maf output from Evolver into a single simulation wide maf. This process is extremely memory intensive with the 120Mb Mammal simulation eventually requiring aprroximately 250Gb of memory.

##Use
1. Write the use section.
2. Follow the use section instructions. ;)
