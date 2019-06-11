## This tutorial covers how to run the [maize python library](https://github.com/orionzhou/maize) on MSI

It is recommended that [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-conda) being installed to handle all dependencies.

Pre-requisites can then be installed using `conda install`:

    conda install biopython numpy pysam pyfaidx

If you plan to run the python scripts on MSI, simply add these lines to your `~/.bashrc` file:

    export maize=/home/springer/zhoux379/git/maize

    for subdir in algorithms apps compara formats projects utils variation test
    do
        if [ -d $maize/$subdir ]; then
            export PATH=$maize/$subdir:$PATH
        fi
    done

Then, sourcing the `~/.bashrc`:

    source ~/.bashrc

Then you should be able to run all the python scripts in the `maize` library:

    (work) zhoux379@ln0005:src $ fasta.py extract -h
    usage: fasta.py extract [-h] [--padding] db loc

    positional arguments:
      db          sequence database (fasta or genome ID)
      loc         location string(s) or BED file(s) (separated by ",")

    optional arguments:
      -h, --help  show this help message and exit
      --padding   padding to size (default: False)

