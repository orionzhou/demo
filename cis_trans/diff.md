* This page shows how to run the cis/trans classification script in a control/treatment experimental setup
  - e.g., whether the gene expression levels before and after stress treatment follows a cis or trans regulatory pattern
  - again we need to obtain raw read counts first, but this time we need the read counts under both the "control" settings and the "treatment" settings
  - [t02_readCount.tsv](t02_readCount.tsv) is a sample read count table with the following columns:
    - `gid`: gene ID
    - `parA_ctrl_rep1`, `parA_ctrl_rep2`, ..., : total read counts in parent A (e.g., B73) replicates under *control* environment
    - `parA_trmt_rep1`, `parA_trmt_rep2`, ..., : total read counts in parent A (e.g., B73) replicates under *treatment* environment
    - `parB_ctrl_rep1`, `parB_ctrl_rep2`, ..., : total read counts in parent B (e.g., Mo17) replicates under *control* environment
    - `parB_trmt_rep1`, `parB_trmt_rep2`, ..., : total read counts in parent B (e.g., Mo17) replicates under *treatment* environment
    - `hybA_ctrl_rep1`, `hybA_ctrl_rep2`, ..., : allele-specific read counts in for allele A (e.g., B73) in hybrid replicates under *control* environment
    - `hybA_trmt_rep1`, `hybA_trmt_rep2`, ..., : allele-specific read counts in for allele A (e.g., B73) in hybrid replicates under *treatment* environment
    - `hybB_ctrl_rep1`, `hybB_ctrl_rep2`, ..., : allele-specific read counts in for allele B (e.g., Mo17) in hybrid replicates under *control* environment
    - `hybB_trmt_rep1`, `hybB_trmt_rep2`, ..., : allele-specific read counts in for allele B (e.g., Mo17) in hybrid replicates under *treatment* environment
  - Parent columns must start with `parA_ctrl_`, `parA_trmt_`, `parB_ctrl_` and `parB_trmt_`, and hybrid columns must start with `hybA_ctrl_`, `hybB_ctrl_`, `hybB_ctrl_` and `hybB_trmt_`
  - Both parents and hybrid can each take one or more replicates
  
* Again, it is recommended to estimate the sample-wise `size factor` and the gene-wise `dispersion` first
  - [t02_sizeFactor.tsv](t02_sizeFactor.tsv) containing `sizeFactors` estimated by DESeq2
  - [t02_dispersion.tsv](t02_dispersion.tsv) containing `dispersions` estimated by DESeq2
  - These two steps are not required, users can set all sample sizeFactors to 1 and gene dispersions to 0 to get the script running

* Run the script in "differential" mode
  ```
  $ ./cis_trans.R --mode diff t01_readCount.tsv t01_sizeFactor.tsv t01_dispersion.tsv t01_result.tsv
  ```
  - example output file: [t02_result.tsv](t02_result.tsv)
    - `gid`: gene ID
    - `prop.p`, `prop.h`: proportion of allele1 expression *change* in the parent and hybrid
    - `reg`: inferred regulatory pattern, one of `conserved`, `cis`, `trans`, `cis+trans` or `unexpected`
  - detailed script usage
  ```
 $ ./cis_trans.R -h
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