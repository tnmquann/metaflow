:description: Quickstart guide to install and run MicrocosM.

Quick Start
============

Install MicrocosM
---------------

Nextflow
~~~~~~~~~~~~

1. Check Java availability

.. code:: bash
   
   java -version

If missing, install Java with SDKMAN:

   a. Install SDKMAN:

   .. code:: bash
   
      curl -s https://get.sdkman.io | bash

   b. Open a new terminal.

   c. Install Java:

   .. code:: bash
   
      sdk install java 17.0.10-tem

   d. Confirm that Java is installed correctly:

   .. code:: bash
   
      java -version

2. Install Nextflow

   a. Download Nextflow:

   .. code:: bash
   
      curl -s https://get.nextflow.io | bash

   b. Make Nextflow executable from everywhere:

   .. code:: bash

      # Make it executable
      chmod +x nextflow
   
      # Move it into an executable path, for example /usr/local/bin
      sudo mv nextflow /usr/local/bin/
   
   .. note::
      Any path that is included in your ``PATH`` variable will work. You can check your ``PATH`` variable by running ``echo $PATH`` in the terminal. If you want to use your own preferred path, you can add it to your ``PATH`` variable by adding the command ``export PATH="/path/to/your/nextflow:$PATH"`` to your bash configuration file, such as ``~/.bashrc`` or ``~/.zshrc``. 

   c. Confirm that Nextflow is installed correctly:

   .. code:: bash
   
      nextflow info

.. hint::
   Follow the `Nextflow documentation <https://www.nextflow.io/docs/latest/install.html>`_
   for full instructions and troubleshooting.

   Also if you're new to Nextflow, we recommend going through the
   `official tutorial <https://training.nextflow.io/latest/nextflow_run/>`_
   for a solid understanding of the platform and its features.

Conda or Mamba
~~~~~~~~~~~~~~~~

If Conda/Mamba is not yet available, install Mamba via `Miniforge <https://github.com/conda-forge/miniforge>`_ following their `installation guide <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`_.

MicrocosM
~~~~~~~~~~~~~~~~

Install through Nextflow portal:

.. code:: bash
   
   nextflow pull tnmquann/metaflow

   # Check pipeline info
   nextflow info tnmquann/metaflow


Prepare taxonomy databases (minimal pre-built versions)
------------------------------------------------------

sourmash
~~~~~~~~~~

Download pre-built database from `ctbrown’s farm <https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/>`_.

We recommend using the `version built on GTDB-RS226 <https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/>`_.

YATCH
~~~~~~~

Download pre-built database available on `Zenodo <https://zenodo.org/communities/yacht/records?q=&l=list&p=1&s=10&sort=newest>`_.

We recommend using the pretrained database on **GTDB-RS226** (ANI threshold 0.995), which is available `here <https://oucruaap-my.sharepoint.com/:f:/g/personal/quantnm_oucru_org/Eu87eLI5JThKg4fNcJEYXAkBUMOSsJB3_MHtNukqVs4oag?e=4tlG0b>`_.

.. important::
   Choose pre-built databases for **sourmash** and **YATCH** that were built on the same **GTDB** or **GenBank** release.


Prepare samplesheet input
------------------

Create a samplesheet file (CSV format) containing paths to your FASTQ files and sample information.

The file must contain at least 5 comma-separated columns, with the following headers:

``sample_id,group,short_reads_1,short_reads_2,long_reads``



Run MicrocosM (on paired-end short reads)
----------------

Read-based workflow
~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
   
   nextflow run tnmquann/metaflow \
      --input /path/to/your/samples.csv \
      --input_format csv \
      -profile conda \
      --outdir /path/to/output/directory \
      --sourmash_database /path/to/your/sourmash_database
      --yacht_database /path/to/your/yacht_database.json \
      --enable_readbase


Assembly-based workflow
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
   
   nextflow run tnmquann/metaflow \
      --input /path/to/your/samples.csv \
      --input_format csv \
      -profile conda \
      --outdir /path/to/output/directory \
      --sourmash_database /path/to/your/sourmash_database

