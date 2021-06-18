## Visualize synteny blocks among syntelogs among multiple species/genotypes

### Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [mamba](https://github.com/mamba-org/mamba)

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    conda install mamba

### Install required software packages (in a new environment):

    conda create -n synteny
    conda env list
    # you should now see a new environment named "synteny"
    conda activate synteny
    mamba install muscle bioconductor-biostrings r-tidyverse r-glue r-ggtree r-tidytree r-ggpubr

### Prepare input files:
- multi-fasta file containing target sequences
  - sequences of the same lengths are recommended
  - [example sequence fasta](test.fa) containing 4kb sequences around TSS (i.e., TSS +/- 2kb) of maize gene "Zm00001d003460" in 25 maize genotypes (i.e., 25 syntenic orthologs using B73 as reference)
- optionally, provide a BED file containing coordinates for highlight
  - [example BED file](test.bed) containing motif coordinates to mark
  - Motif positions can be obtained using a wrapper script [fimo.py](https://github.com/orionzhou/nf/blob/master/bin/mmm/fimo.py) that runs [FIMO](https://meme-suite.org/meme/tools/fimo) internally:

    fimo.py locate --motif b01_l002_5 /home/springer/zhoux379/stress/data/41_ml/00_nf/03_motif_lists/b01.meme test.fa test.bed

- R data file containing gene structure information for provided orthologs
  - pre-made file available for maize B73v4-anchored ortholog structures: `/home/springer/zhoux379/projects/genome/data2/syntelog/maize.genes.v4.rds`

Finally, run the synteny visualization script:

    ./syntelog.R --gid Zm00001d003460 test.fa test.bed test.pdf

For detailed help:

```bash
(synteny) pzhou@sugon1:syntelog $ ./plot_syntelog.R -h
usage: ./plot_syntelog.R [-h] [--gid GID] [--gene GENE] [--width WIDTH]
                         [--height HEIGHT]
                         f_fa f_bed fp

Make synteny plot using provided fasta sequences and mark features

positional arguments:
  f_fa             multi-fasta sequences
  f_bed            BED file containing feature positions to highlight
  fp               output (PDF) file

optional arguments:
  -h, --help       show this help message and exit
  --gid GID        gene ID used to extract structure info [default: test gene]
  --gene GENE      R data file containing gene structure information for
                   provided orthologs [default: /datalus/weiyu/projects/genome
                   /data2/syntelog/maize.genes.v4.rds]
  --width WIDTH    figure width [default: 7]
  --height HEIGHT  figure height [default: 6]
```