# Batch design primers using Primer3

## Requirements
1. Install [Conda](https://docs.conda.io/en/latest/miniconda.html) (if you already have conda installed, skip to Step 2):

        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh

2. Install required packages

        # create a conda environment
        conda create -n primer3

        # switch to the newly created enrironment
        conda activate primer3

        # install required packages
        conda install primer3 primer3-py pyyaml pandas

3. Get the script and primer3 configuration file:

        git clone git@github.com:orionzhou/demo.git
        cd demo/primer3

  * The primer3 global configuration can be found [primer3.yml](primer3.yml)

## Align RNA-Seq reads to the reference

