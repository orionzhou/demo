# Classify cis/trans inheritance pattern using inbred/hybrid expression data

* Install R packages: `tidyverse`, `bbmle`, `multidplyr`
	
	conda install r-tidyverse r-bbmle r-multidplyr

* To run the cis/trans classification script, the first thing we need to do is to obtain raw read counts for our samples.
  - Specifically, we need to obtain total read counts for the parents (e.g., B73 and Mo17 inbred lines) as well as allele-specific read counts for both allels in the hybrid
  - Although multiple replicates are not required to run the script, >=2 replicates are highly recommended since it significantly increases statistical confidence
  - [t01_readCount.tsv](t01_readCount.tsv) is a sample read count table with the following columns:
    - `gid`: gene ID
    - `parA_rep1`, `parA_rep2`, ..., : total read counts in parent A (e.g., B73) replicates
    - `parB_rep1`, `parB_rep2`, ..., : total read counts in parent B (e.g., Mo17) replicates
    - `hybA_rep1`, `hybA_rep2`, ..., : allele-specific read counts in for allele A (e.g., B73) in hybrid replicates
    - `hybB_rep1`, `hybB_rep2`, ..., : allele-specific read counts in for allele B (e.g., Mo17) in hybrid replicates
  - Parent columns must start with `parA_` and `parB_`, and hybrid columns must start with `hybA_` and `hybB_`
  - Both parents and hybrid can each take one or more replicates
  
* Before running the script it is also recommended to estimate two parameters using DESeq2, the sample-wise `size factor` as well as the gene-wise `dispersion`
  - check this excellent walkthrough and explanation on [what these two parameters mean and how to estimate them with DESeq2](https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)
  - example sample-wise size factor table [t01_sizeFactor.tsv](t01_sizeFactor.tsv) containing `sizeFactors` estimated by DESeq2
  - example gene-wise dispersion table [t01_dispersion.tsv](t01_dispersion.tsv) containing `dispersions` estimated by DESeq2
  - These two steps are not required, users can set all sample sizeFactors to 1 and gene dispersions to 0 to get the script running

* Run the script
  
  $ ./cis_trans.R t01_readCount.tsv t01_sizeFactor.tsv t01_dispersion.tsv t01_result.tsv

  - example output file: [t01_result.tsv](t01_result.tsv)
    - `gid`: gene ID
    - `prop.p`, `prop.h`: proportion of allele1 expression in the parent and hybrid
    - `reg`: inferred regulatory pattern, one of `conserved`, `cis`, `trans`, `cis+trans` or `unexpected`

```
 $ ./cis_trans.R -h
NULL
Warning message:
package ‘argparse’ was built under R version 4.0.3
usage: ./cis_trans.R [-h] [--mode MODE] [--min_rc MIN_RC] [--n_cpu N_CPU]
                     f_rc f_sf f_dsp fo

Classify cis/trans inheritance pattern using inbred/hybrid RNA-Seq read counts

positional arguments:
  f_rc             read count table
  f_sf             sample-wise size factor table
  f_dsp            gene-wise dispersion table
  fo               output file

optional arguments:
  -h, --help       show this help message and exit
  --mode MODE      cis/trans test mode, "basic" for steady-state cis/trans
                   test and "diff" for control/treatment differential test
                   [default: basic]
  --min_rc MIN_RC  minimum read counts to filter low-expressed genes [default:
                   10]
  --n_cpu N_CPU    number of CPUs / threads to use for parallel processing
                   (for spped up if you have many genes) [default: 1]
```