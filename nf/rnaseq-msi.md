## Nextflow RNA-Seq pipeline

### Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [mamba](https://github.com/mamba-org/mamba)

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    conda install mamba

### Intall [Nextflow](https://github.com/nextflow-io/nextflow) (in a new conda environment)

    conda create -n nf python=3
    conda activate nf
    conda install nextflow
    conda list

### Clone the repo to local:

    cd /home/springer/zhoux379/git
    git clone https://github.com/orionzhou/nf.git

(or if you cloned the repo a while ago) update the local repo (if there are new changes):

    cd /home/springer/zhoux379/git/nf
    git pull

### Create a conda environment named `rnaseq` to run the pipeline:
This can be done either [using an existing environment yml file](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file):

    conda env create -n rnaseq -f $NXF_HOME/configs/environments/rnaseq.yml
    conda env list
    # there should now be a new environment named "rnaseq"

or create an empty environment and then manual install the packages:

    conda create -n rnaseq
    conda env list
    # there should now be a new environment named "rnaseq"
    conda activate rnaseq
    mamba install parallel sra-tools pigz fastqc trim-galore bwa hisat2 star picard samtools bcftools bedtools bamtools pysam sambamba preseq alfred bioawk biopython deeptools kallisto minimap2 numpy pandas perl plotly pyfaidx pybigwig ucsc-bedgraphtobigwig subread r-base r-tidyverse r-glue r-argparse r-purrr r-readr r-readxl rseqc bioconductor-dupradar bioconductor-edger bioconductor-deseq2 salmon stringtie multiqc

### Post set-up
Add these environmental variables (with necessary modification) to your `~/.bashrc` or `~/.bash_profile`

    export NXF_HOME=/home/springer/zhoux379/git/nf
    export NXF_CACHE=/scratch.global/zhoux379/nf
    export NXF_CONDA_CACHEDIR=/home/springer/zhoux379/software/miniconda3/envs
    export NXF_WORK=$NXF_CACHE/work
    export NXF_TEMP=$NXF_CACHE/tmp
    export NXF_OPTS='-Xms1g -Xmx10g'

Log out and log in again (or run `source ~/.bashrc`) to make these variables into effect

Copy the test pipeline (`$NXF_HOME/test/rnaseq`) to wherever you'd like to take a try:

    cp -rf $NXF_HOME/test/rnaseq /home/springer/zhoux379/rnaseq_test
    cd ~/rnaseq_test
    ls
    # genomes.yml nextflow.config  design.tsv  design.xlsx

Make necessary changes to:
- `design.tsv` which contains paths to your fastq sequences
- `nextflow.config`: in particular, `launchDir=XXX` and `email_on_fail=123@gmail.com`
- as well as the genome index configuration file `genomes.yml`

Finally, (inside this directory) we are set to run the test pipeline using existing genome database on MSI-mesabi queue:

    nextflow run $NXF_HOME/rnaseq -params-file genomes.yml -profile conda
