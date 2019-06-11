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
  * The actual script can be found [run_primer3.py](run_primer3.py)

## Running the script

* Input file:
  * the input needs to a tab-separated text file with no headers and 4 columns:
    * column 1: an arbitrary ID for the sequence, need to be unique
    * columns 2 and 3: **relative** location (start, end) of the target segment in the sequence
    * column 4: input sequence
  * A sample input file can be found here: [sample.tsv](sample.tsv)
* run the script

        # show help message
        ./run_primer3.py -h
        ./run_primer3.py ufmu -h

        # run the script with sample input
        ./run_primer3.py ufmu sample.tsv > sample_out.tsv

* Output
  * the output file will have the following columns:
    * `pid`, `start`, `end`: the unique sequence ID, target start and end positions as in the input
    * `n_left`, `n_right`, `n_pair`: number of left/right primers and primer pairs returned by Primer3
    * `left.start`, `left.size`, `left.seq`, `right.start`, `right.size`, `right.seq`: the start position, size and sequence of the left/right primer
  * a sample output file can be found here: [sample_out.tsv](sample_out.tsv)


