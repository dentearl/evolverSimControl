#evolverSimControl

[Dent Earl](https://github.com/dentearl/)
2009 - 2011

A [jobTree](https://github.com/benedictpaten/jobTree/) based wrapper for the [Evolver](http://www.drive5.com/evolver/) genome evolution simulation tool suite.

##Dependencies
* sonLib: https://github.com/benedictpaten/sonLib/
* jobTree: https://github.com/benedictpaten/jobTree/
* Evolver: http://www.drive5.com/evolver/ Specifically, make sure that the Evolver tools are on your <code>PATH</code> environmental variable and that their names are preceeded with <code>evolver_</code>. Specifically all of the following list of files need to be on your <code>PATH</code>. Depending on how you add these things to your <code>PATH</code>, you may need to edit some of the evolver python scripts to ensure that they can import the evolver <code>gff.py</code> library.
    * <code>evolver_codon_report.pl</code>
    * <code>evolver_cvt</code>
    * <code>evolver_drawrev</code>
    * <code>evolver_evo</code>
    * <code>evolver_evostats_report.py</code>
    * <code>evolver_gff_cdsutr2exons.py</code>
    * <code>evolver_gff_exons2introns.py</code>
    * <code>evolver_gff_featurestats2.py</code>
    * <code>evolver_gff_featurestats2.sh</code>
    * <code>evolver_handle_mobiles.pl</code>
    * <code>evolver_gene_deactivate.sh</code>
    * <code>evolver_merge_evostats.py</code>
    * <code>evolver_mobile_report.pl</code>
    * <code>evolver_transalign</code>
    * <code>evolver_trf2gff.py</code>
* mafJoin: https://github.com/dentearl/mafTools Not necessary for simple simulations, mafJoin (part of mafTools) is only needed if you wish to create an maf alignment of all sequences following a simulation.

##Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.
4. Edit your <code>PYTHONPATH</code> variable to contain the parent directory of the <code>evolverSimControl/</code> directory.
5. Type <code>make test</code>.

##Use

##Example
