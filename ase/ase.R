### these packages need to be installed first
require(VGAM)
require(gamlss)
require(BiocParallel)
require(tidyverse)

### Parameters
# sample input file (demo_parents.tsv) has 6 columns:
#   gene ID, Genotype (B73 or Mo17), Replicate (at least one),
#   raw read count, (edgeR-normalized) CPM
input_parent = './demo_parents.tsv'

# sample input file (demo_hybrid.tsv) also has 6 columns:
#   gene ID, Genotype (BxM), Replicate (at least one),
#   raw allele specific read counts supporting the
#   reference allele (nref), alternate allele (nalt) and conflict reads (ncft)
#   in the F1 hybrid
input_hybrid = './demo_hybrid.tsv'

## filters to remove bad genes
# minimum number of allele specific reads
min_asr = 10
# maximum proportion of conflicting reads
max_cft = 0.05

# number of CPUs / threads to use for parallel processing (if you have many genes)
n_cpu = 1


### load functions from "ase_functions.R"
source("ase_functions.R")

### run the function to classify cis/trans patterns
ase = classify_cis_trans(input_parent, input_hybrid, min_asr, max_cft, n_cpu)

# save to output file
fo = './demo_out.tsv'
write_tsv(ase, fo)
