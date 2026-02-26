:description: Home

microcosM
=============

.. rst-class:: lead

   Optimised, Comprehensive Metagenomics

----

microcosM is a modular, robust pipeline for streamlining complex metagenomic data analysis workflows.

microcosM provides solutions to both rapid, read-level taxonomic classification and high-resolution genome assembly. With a modular design and a rich suite of updated, carefully curated tools, data analysis workflows can be easily customised in a trackable manner to meet diverse research questions and resource settings.


Key Features
----------------------

Read-based Analysis
~~~~~~~~~~~~~~~~~~~~~~~

- Taxonomic assignment for each sequencing read is performed using an optimized combination of sourmash_ and YACHT_. This approach ensures fast, scalable processing, even for large datasets. Reference databases can be customized to fit your specific research needs.
- Detection of antimicrobial resistance genes (ARGs) is integrated via rgi_bwt_, providing actionable insights into resistance profiles.

.. _sourmash: https://github.com/sourmash-bio/sourmash
.. _YACHT: https://github.com/KoslickiLab/YACHT
.. _rgi_bwt: https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst

Metagenome Assembly & Binning
~~~~~~~~~~~~~~~~~~~~~~~

- The pipeline assembles contigs and performs binning to recover Metagenome Assembled Genomes (MAGs), enabling in silico characterization of microbial communities.


.. toctree::
    :caption: Getting started
    :hidden:

    quickstart
    installation/index
    
