#evolverSimControl

[Dent Earl](https://github.com/dentearl/)
2009 - 2011

A [jobTree](https://github.com/benedictpaten/jobTree/) based wrapper for the [Evolver](http://www.drive5.com/evolver/) genome evolution simulation tool suite.

##Dependencies
* jobTree: https://github.com/benedictpaten/jobTree/
* evolver: http://www.drive5.com/evolver/
* mafJoin: https://github.com/dentearl/mafTools Not necessary for simulation, only needed if you wish to create an alignment of all sequences following a simulation.

##Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.
4. Make sure that the _Evolver_ tools are on your <code>PATH</code> and that their names are preceeded with <code>evolver_</code>. Specifically the following need to be on your <code>PATH</code>:
        * <code>evolver_cdsutr2exons.py</code>
        * <code>evolver_codon_report.py</code>
        * <code>evolver_cvt</code>
        * <code>evolver_drawrev</code>
        * <code>evolver_evo</code>
        * <code>evolver_evostats_report.py</code>
        * <code>evolver_exons2introns.py</code>
        * <code>evolver_gff_featuresats2.py</code>
        * <code>evolver_handle_mobiles.pl</code>
        * <code>evolver_handle_gene_deactivate.sh</code>
        * <code>evolver_merge_evostats.py</code>
        * <code>evolver_mobile_report.py</code>
        * <code>evolver_transalign</code>
        * <code>evolver_trf2gff.py</code>
5. Depending on how you add these things to your <code>PATH</code> need to edit some of the evolver python scripts to ensure that they can import the evolver <code>gff.py</code> library.

##Use
