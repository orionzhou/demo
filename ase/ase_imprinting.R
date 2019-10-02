### these packages need to be installed first
require(VGAM)
#require(BiocParallel)
require(tidyverse)

### Parameters
# sample input file (demo_hybrid.tsv) also has 5 columns:
#   gene ID, Genotype (BxM), Replicate (at least one),
#   raw allele specific read counts supporting the allele1 and allele2
#   in the F1 hybrid
input_hybrid = './demo_imprinting.tsv'

## filters to remove bad genes
# minimum number of allele specific reads
min_asr = 10

# read ASE read counts
ti = read_tsv(input_hybrid)

# prepare data
ta = ti %>% rename(n1 = allele1, n2 = allele2) %>%
    group_by(gid) %>%
    summarise(nl1 = list(n1), nl2 = list(n2), n1=sum(n1), n2=sum(n2)) %>%
    ungroup() %>%
    filter(n1 + n2 >= min_asr)

# beta-binomail likelihood function
fit_bb <- function(nl1, nl2) {
    #{{{
    # a priori frequency
    prob = 2/3
    # total number of reads
    x = nl1; size = nl1 + nl2
    # log likelihood function
    LL <- function(prob, rho) {
        if(prob > 0 & prob < 1 & rho >= 0 & rho < 1)
            -sum(VGAM::dbetabinom(x, size, prob = prob, rho = rho, log = TRUE))
        else 100
    }
    # number of replicates
    nrep = length(nl1)
    # fit the BB model centering around 0.67 (degree of freedome = 1)
    fit1 = mle(LL, start = list(prob = prob, rho = 0),
          method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5),
          fixed = list(prob = prob),
          nobs = nrep)
    # fit a free BB model with 2 degrees of freedom
    fit2 = mle(LL, start = list(prob = prob, rho = 0),
          method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5),
          nobs = nrep)
    bic = BIC(fit1, fit2)
    bic1 = bic$BIC[1]; bic2 = bic$BIC[2]
    tibble('bic1'=bic1, 'bic2'=bic2)
    #}}}
}

# bic1 and bic2 are the -log(likelihood) of model fitting observed data
# smaller bic means more likely model
tb = ta %>% mutate(res = map2(nl1, nl2, fit_bb)) %>%
    select(gid, n1, n2, res) %>% unnest() %>%
    mutate(imprint = ifelse(bic1 < bic2, 'no_imprint', 'imprint'))

# currently the script calls an impriting whenever model1 outperforms model2 (bic1 - bic2 > 0)
# can also require bic1 - bic2 >= 2 to make a more stringent call

