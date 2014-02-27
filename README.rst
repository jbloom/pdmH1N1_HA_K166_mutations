----------------------------------------------------------
Mutations at site 166 in HA of human pandemic H1N1
----------------------------------------------------------

Analysis of mutations at site 166 (H3 numbering, see `HA numbering`_ below) of the HA from human pandemic H1N1. 

The input files, analysis scripts, and results can be downloaded `on GitHub`_.

This analysis was performed by `Jesse Bloom`_.

.. contents:: Contents
   :depth: 2

Methods summary
----------------
The occurrence of different amino-acid identities at HA residue 166 (H3 numbering) was analyzed by downloading all full-length human pandemic H1N1 sequences present in the `Influenza Virus Resource`_ [Bao2008]_ as of Feb-23-2014. After purging sequences that were less than full length, contained ambiguous nucleotide identities, lacked full (year, month, day) isolation dates, or were otherwise anomalous, the sequences were aligned. Each calendar year was broken into four equal partitions beginning with January 1, the frequencies of different amino acids at residue 166 for each partition was calculated and plotted. For construction of phylogenetic trees, the sequence set was randomly subsampled to 10 sequences per quarter-year partition. 
`BEAST`_ [Drummond2013]_ was then used to sample from the posterior distribution of phylogenetic trees with reconstructed sequences at the nodes, after date stamping the sequences, using a `JTT`_ [Jones1992]_ with a single rate category with an exponential prior, a strict molecular clock, and relatively uninformative coalescent-based prior over the tree. The figure shows a maximum clade credibility summary of the posterior distribution with branches colored according to the reconstructed amino-acid identity with the highest posterior probability at their descendent nodes. The tree was visually rendered using `FigTree`_. The input data and computer code used for this analysis can be found `on GitHub`_ at https://github.com/jbloom/pdmH1N1_HA_K166_mutations.

HA numbering
-------------
Various numbering schemes are in use for influenza HA (see `HA_numbering`_). The site examined here is:

* Residue 166 in the H3 numbering scheme (for example, as in PDB `4HMG`_).

* Residue 169 in the H1 numbering scheme (for example, as in PDB `4JTV`_).

* Residue 180 in sequential 1, 2, ... numbering of typical pandemic H1N1 HAs beginning with 1 a the N-terminal methionine.

Software used
---------------
* `Python`_ : version 2.7.5

* `EMBOSS needle`_ : from EMBOSS 6.3.1, for sequence alignments.

* `mapmuts`_: version 1.0, for sequence manipulation utilities.

* `matplotlib`_ : version 1.3.1, for plotting.

* `ImageMagick convert`_ : version 6.8.3-3, for converting PDF images to JPGs.

* `BEAST`_ : version 1.8.0, for inferring phylogenetic tree. Executables used from the package are ``beast`` and ``treeannotator``.

* `BEAGLE`_ : version 2.1.2, for use with `BEAST`_.

* `Tracer`_ : version 1.6, for analyzing `BEAST`_ output.

* `FigTree`_ : version 1.4, for creating the tree image.

Input files
-------------
Here are the input sequence data files. All other files are analysis scripts or output files that form part of the `Analysis`_ described below. Here are the input files:

* *HAseqs.fasta* : File containing all full-length protein-coding human pandemic H1N1 HA sequences as downloaded from the `Influenza Virus Resource`_ on February-23-2014.

* *California2009-HA.fasta* : File containing protein-coding HA from *A/California/4/2009*, which is used as the reference sequence for alignments.

* *AnomalousSequences.txt* : File listing sequences that appear to be anomalous due to much higher than expected divergence. These are probably sequences that were incorrectly classified as human pandemic H1N1 in the `Influenza Virus Resource`_. These sequences are excluded from the analyses.

Analysis
-----------
Here are the steps in the analysis:

Sequence alignment, sequence selection, mutation counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The first step is performed by the `Python`_ script *parseseqs.py*, and can be run with::

    python parseseqs.py

The steps are:

1) Align the sequences. All sequences in *HAseqs.fasta* are pairwise aligned with the reference sequence in *California2009-HA.fasta* after:

    a) Purging any sequences with ambiguous nucleotides, and any sequences that are specified as anomalies in *AnomalousSequences.txt*.

    b) Purging any sequences that do not exactly match the reference sequence in length.

    c) Purging any sequences that do not have a full date specified (year, month, and day).

  The resulting aligned and translated protein sequences are written to the created file *all_alignedproteins.fasta*.

2) Get a subset of sequences for building the phylogenetic tree. Each year is broken into four partitions spanning months 1-3, 4-6, 7-9, and 10-12. For each partition, up to 10 sequences per year are randomly selected to be retained. These sequences are written to *selected_alignedproteins.fasta*. The reason for only keeping a subset of the sequence in this way is to have a sequence set small enough to build a phylogenetic tree without undue computational burden.

3) Count the occurrences of different amino-acid identities at site 166 (H3 numbering, see `HA numbering`_) for each of the date partitions. The results are written to the file *partitioncounts.txt*. Here is the raw data in that file:

    .. include:: partitioncounts.txt
       :literal:

Plotting amino-acid fractions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The `Python`_ script *makeplot.py* is used to make a PDF plot (with `matplotlib`_) showing the frequencies of different amino acids at site 166 as a function of time. Run the script with::

    python makeplot.py

The resulting plot divides each year into four partitions (months 1-3, 4-6, 7-9, 10-12) as in *partitioncounts.txt*. The frequencies of different amino acids at site 166 is then shown for each of these partitions.

The created plot is *aafracs.pdf*. A JPG version (lower quality, created with `ImageMagick convert`_ from the PDF) is also created as *aafracs.jpg*. Here is that plot:

  .. figure:: aafracs.jpg
     :alt: aafracs.jpg
     :width: 40%
     :align: center

     The plot *aafracs.jpg* showing the frequencies of different amino acids at site 166. A higher quality image of this same plot is in *aafracs.pdf*.

This plot is a graphical display of the data in *partitioncounts.txt*.


Construction of phylogenetic tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A phylogenetic tree was constructed using the subset of sequences in *selected_alignedproteins.fasta*. This subset only contains 10 sequences per year partition (4 partitions per year) -- the reason for only using some sequences is to make the phylogenetic tree construction tractable.

First, the ``beauti`` program of the `BEAST`_ package was used to construct an input XML file from *selected_alignedproteins.fasta*. The tips were stamped with their isolation dates, a `JTT`_ substitution model with a single rate category, and a strict molecular clock was assumed. The tree prior was a relatively uninformative coalescent-based prior and an exponential prior was used over the rate. Ancestral states were reconstructed at all nodes. The MCMC was run for 10 million steps. Full details are available in the created XML file *selected_alignedproteins.xml*.

`BEAST`_ was then run using `BEAGLE`_ with the command::

    ~/BEASTv1.8.0/bin/beast -overwrite -beagle_SSE -seed 1 -threads 3 selected_alignedproteins.xml

    java -classpath ~/BEASTv1.8.0/lib/beast.jar dr.app.beast.BeastMain -threads 6 -overwrite -beagle_SSE selected_alignedproteins.xml

This created the following files::

    selected_alignedproteins.log
    selected_alignedproteins.trees

The key data is in the latter file (*selected_alignedproteins.trees*); however, this file is very large and so is not included in this repository on `on GitHub`_ (it can be regenerated using the commands above).

`Tracer`_ was used to analyze the *selected_alignedproteins.log* file to check for MCMC convergence. If the first 10% of steps (first 1 million steps) are excluded, the remaining 90% (9 million steps) appear to be well converged, and have good effective sample sizes. These results indicate that the MCMC was probably run for a sufficiently large number of steps.

An annotated tree was then constructed using ``logcombiner`` and ``treeannotator`` from the `BEAST`_ package. This tree is created automatically by the `Python`_ script *analyze_tree.py* with the command::

    python analyze_tree.py

This script creates several output files not included in the package `on GitHub`_ (again, they can be regenerated) -- the script also contains hardcoded paths to `BEAST`_ executables that you may need to modify for your computer. The result of this script that is retained is the file *annotated_maxcredtrees.trees*.

The *annotated_maxcredtrees.trees* file was then opened with `FigTree`_ and manually re-formatted for appealing visual appearance. The formatted image was saved as *annotated_maxcredtree.pdf*. A JPG version of this file (*annotated_maxcredtree.jpg*) was then created with `ImageMagick convert`_ using::

    convert -density 300 annotated_maxcredtree.pdf annotated_maxcredtree.jpg

Here is that plot:

  .. figure:: annotated_maxcredtree.jpg
     :alt: annotated_maxcredtree.jpg
     :width: 80%
     :align: center

     The plot *annotated_maxcredtree.jpg* showing the frequencies of different amino acids at site 166. A higher quality image of this same plot is in *annotated_maxcredtree.pdf*. In this plot, the residue is labeled as 180 since `BEAST`_ uses consecutive numbering rather than the H3 numbering scheme.



.. _`Influenza Virus Resource`: https://www.ncbi.nlm.nih.gov/genomes/FLU/FLU.html
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Python`: http://www.python.org/
.. _`EMBOSS needle`: http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
.. _`HA_numbering`: https://github.com/jbloom/HA_numbering
.. _`4HMG`: http://www.rcsb.org/pdb/explore.do?structureId=4HMG
.. _`4JTV`: http://www.rcsb.org/pdb/explore.do?structureId=4JTV
.. _`mapmuts`: http://jbloom.github.io/mapmuts/
.. _`on GitHub`: https://github.com/jbloom/pdmH1N1_HA_K166_mutations
.. _`matplotlib`: http://matplotlib.org/
.. _`ImageMagick convert`: http://www.imagemagick.org/script/index.php
.. _`BEAST`: http://beast.bio.ed.ac.uk/Main_Page
.. _`JTT`: http://www.ncbi.nlm.nih.gov/pubmed/1633570
.. _`BEAGLE`: http://beast.bio.ed.ac.uk/BEAGLE
.. _`Tracer`: http://beast.bio.ed.ac.uk/Tracer
.. _`FigTree`: http://tree.bio.ed.ac.uk/software/figtree/
.. [Drummond2013] Drummond AJ, Suchard MA, Xie D, and Rambaut A. Bayesian phylogenetics with BEAUti and BEAST 1.7. Mol Biol Evol. 29:1969-1973 (2012)
.. [Bao2008] Bao Y, Dernovoy D, Kiryutin B, Zaslavsky L, Tatusova T, Ostell J, and Lipman D. The influenza virus resource at the National Centery for Biotechnology Information. J Virol. 82:596-601 (2008).
.. [Jones1992] Jones DT, Taylor WR, and Thornton JM. The rapid generation of mutation data matrices from protein sequences. Comput Appl Biosci. 8:275-282 (1992)
