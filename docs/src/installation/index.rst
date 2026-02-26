Installation
=============

Requirements
----------------

Software
~~~~~~~~~~~~~~

* Nextflow (installation guide `here <https://www.nextflow.io/docs/latest/install.html>`_)
* Conda or Mamba (installation guide for Mamba (recommended) `here <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`_)
* Optional: A container runtime engine like `Docker <https://www.docker.com/>`_ or `Apptainer <https://apptainer.org/docs/admin/main/installation.html>`_.

Hardware
~~~~~~~~~~~~

* A POSIX-compatible operating system (e.g., Linux, macOS) or Windows with `WSL <https://learn.microsoft.com/en-us/windows/wsl/>`_. We strongly recommend using a Unix-based OS for best compatibility with all tools.
* At least **32 GB of RAM** for **MEGAHIT**, or **128 GB** if you plan to use **metaSPAdes** for assembly.
* At least **256 GB of disk space**.

.. admonition:: ℹ️ Storage Estimates (based on shallow shotgun metagenomics)

   - Nextflow installation: ~150 MB
   - Databases: ~120 GB
   - Environment folder: ~20 GB
   - Intermediate/cache files: 20–50 GB (varies with sample size and quality)
   - Output per sample: ~9 GB

* *Optional*: Some steps in assembly-based analysis (for example, **COMEBin** in the binning step) will run faster if your system is equipped with a GPU.

Install MicrocosM
-----------------

via Nextflow portal (recommended)
~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
   
   nextflow pull tnmquann/metaflow

   # Check pipeline info
   nextflow info tnmquann/metaflow

by cloning the repository (dev versions available)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
   
   git clone https://github.com/tnmquann/metaflow.git
   cd metaflow
   nextflow run main.nf 

Database Preparation
-----------------

* Minimal
* Extended
