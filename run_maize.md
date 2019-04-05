## This tutorial covers how to run the [maize python library](https://github.com/orionzhou/maize) on MSI

It is recommended that [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-conda) is installed to handle all dependencies.

Pre-requisites can then be installed using `conda install`:

    conda install biopython numpy pysam pyfaidx

If you plan to run the python scripts on MSI, you simply need to add these lines to your `~/.bashrc` file:

    export maize=/home/springer/zhoux379/git/maize

    for subdir in algorithms apps compara formats projects utils variation test
    do
        if [ -d $maize/$pname ]; then
            export PATH=$maize/$subdir:$PATH
        fi
    done

Then, sourcing the `~/.bashrc`:

    source ~/.bashrc

Then you should be able to run all the python scripts in the `maize` library:

    fasta.py -h

