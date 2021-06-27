## Visualize synteny blocks among syntelogs among multiple species/genotypes

### Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [mamba](https://github.com/mamba-org/mamba)

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    conda install mamba

### Install required software packages (in a new environment):

    conda create -n synteny
    conda env list
    # there should now be a new environment named "synteny"
    conda activate synteny
    mamba install mafft bioconductor-biostrings r-tidyverse r-glue bioconductor-ggtree r-tidytree r-ape r-ggpubr

### Prepare input files:
- multi-fasta file containing target sequences
  - sequences of the same lengths are recommended but not required
  - [example sequence fasta](test2.fa) containing 4kb sequences around TSS (i.e., TSS +/- 2kb) of maize gene "Zm00001d003460" in 25 maize genotypes (i.e., 25 syntenic orthologs using B73 as reference)
  - use [this script](prep_seq.md) to extract syntelog sequences from multiple genome assemblies
- optionally, provide a BED file containing coordinates for highlight
  - [example BED file](test2.bed) containing motif coordinates to mark
  - Motif positions can be obtained using a wrapper script [fimo.py](https://github.com/orionzhou/nf/blob/master/bin/mmm/fimo.py) that runs [FIMO](https://meme-suite.org/meme/tools/fimo) internally:
    ```
    fimo.py locate --motif b01_l002_5 /home/springer/zhoux379/stress/data/41_ml/00_nf/03_motif_lists/b01.meme test.fa test.bed
    ```
- R data file containing gene structure information for provided orthologs
  - pre-made file available for:
    - maize B73v4-anchored ortholog structure in 30 maize genotypes: `/home/springer/zhoux379/projects/genome/data2/syntelog/xref.maize.v4.rds`
    - maize B73v5-anchored ortholog structure in 30 maize genotypes: `/home/springer/zhoux379/projects/genome/data2/syntelog/xref.maize.v5.rds`
    - wheat CS-D-anchored ortholog structure in 30 maize genotypes: `/datalus/weiyu/projects/genome/data2/syntelog/xref.wheat.csd.rds`

Finally, run the synteny visualization script:

    ./plot_syntelog.R --gid Zm00001d003460 --bed test2.bed test2.fa test2.pdf

For detailed help:

```bash
(synteny) pzhou@sugon1:syntelog $ ./plot_syntelog.R -h
usage: ./plot_syntelog.R [-h] [--gid GID] [--opt OPT] [--bed BED]
                         [--gene GENE] [--width WIDTH] [--height HEIGHT]
                         fi fo

Make synteny plot using provided fasta sequences and mark features

positional arguments:
  fi               multi-fasta sequences
  fo               output (PDF) file

optional arguments:
  -h, --help       show this help message and exit
  --gid GID        gene ID used to extract structure info [default: test gene]
  --opt OPT        range options [default: genomic+2k]
  --bed BED        BED file containing feature positions to highlight
  --gene GENE      R data file containing gene structure information for
                   provided orthologs [default: /home/springer/zhoux379/projec
                   ts/genome/data2/syntelog/xref.maize.v4.rds]
  --width WIDTH    figure width [default: 7]
  --height HEIGHT  figure height [default: 6]
```
