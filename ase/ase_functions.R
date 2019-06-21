require(VGAM)
require(gamlss)
require(BiocParallel)
require(tidyverse)

get_bic <- function(i, tt) {
    #{{{
    prob = tt$B73[i] / (tt$B73[i] + tt$Mo17[i])
    x = unlist(tt$refl[i]); size = unlist(tt$refl[i]) + unlist(tt$altl[i])
    prob.h = sum(x) / sum(size)
    bicc = NA; bict = NA; bicct = NA; reg = NA
    if(prob == 0)
        reg = ifelse(prob.h < 0.05, 'cis',
              ifelse(prob.h > 0.45 & prob.h < 0.55, 'trans', 'cis+trans'))
    else if(prob == 1)
        reg = ifelse(prob.h > 0.95, 'cis',
              ifelse(prob.h > 0.45 & prob.h < 0.55, 'trans', 'cis+trans'))
    else {
        LL <- function(prob, rho) {
            if(prob > 0 & prob < 1 & rho >= 0 & rho < 1)
                -sum(VGAM::dbetabinom(x, size, prob = prob, rho = rho, log = TRUE))
            else 100
        }
        fitc = mle(LL, start = list(prob = prob, rho = 0),
              method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5),
              fixed = list(prob = prob),
              nobs = length(x))
        fitt = mle(LL, start = list(prob = prob, rho = 0),
              method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5),
              fixed = list(prob = 0.5),
              nobs = length(x))
        fit = mle(LL, start = list(prob = prob, rho = 0),
              method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5),
              nobs = length(x))
        #coef(fitc)
        bic = BIC(fitc, fitt, fit)
        bicc = bic$BIC[1]; bict = bic$BIC[2]; bicct = bic$BIC[3]
        tb = as_tibble(bic) %>%
            add_column(reg = c('cis', 'trans','cis+trans')) %>% arrange(BIC)
        reg = tb$reg[1]
    }
    c('cis'=bicc, 'trans'=bict, 'cis+trans'=bicct, 'reg'=reg)
    #}}}
}

classify_cis_trans <- function(input_parent, input_hybrid, min_asr=10, max_cft=.05, n_cpu=1) {
    #{{{ MLE approach to estimate cis/trans
    ti1 = read_tsv(input_parent)
    ti2 = read_tsv(input_hybrid)
    # you can always check the result by typing the name of the data frame / tibble such as "ti1"
#
    bpparam <- MulticoreParam(n_cpu)
    ta = ti2 %>% rename(n0 = nref, n1 = nalt) %>%
        group_by(gid) %>%
        summarise(ref = sum(n0), alt = sum(n1), cft = sum(ncft),
               refl = list(n0), altl = list(n1), cftl = list(ncft)) %>%
        mutate(ntc = ref + alt + cft, nt = ref + alt, pcft = cft / ntc) %>%
        filter(pcft <= max_cft, nt >= min_asr) %>%
        ungroup() %>% dplyr::select(-ntc, -nt, -pcft)
    tt = ti1 %>% group_by(gid, Genotype) %>%
        summarise(CPM=mean(CPM)) %>%
        spread(Genotype, CPM) %>%
        inner_join(ta, by = 'gid')
#
    y = bplapply(1:nrow(tt), get_bic, tt, BPPARAM = bpparam)
    tb = do.call(rbind.data.frame, y) %>% as_tibble()
    colnames(tb) = c("bic.c", "bic.t", "bic.ct", "reg")
    tb = tb %>% mutate(bic.c = as.numeric(bic.c),
                       bic.t = as.numeric(bic.t),
                       bic.ct = as.numeric(bic.ct))
    stopifnot(nrow(tb) == nrow(tt))
    ttb = tt %>% bind_cols(tb) %>%
        dplyr::select(-refl, -altl, -cftl) %>%
        as_tibble()
#
    reg_levels = c("cis only", "trans only", "cis+trans", "cis B73", "cis Mo17", "trans B73", "trans Mo17")
    t_ase = ttb %>%
        mutate(prop.p = B73/(B73+Mo17), prop.h = ref/(ref+alt)) %>%
        mutate(reg = ifelse(reg == 'cis', 'cis only',
                     ifelse(reg == 'trans', 'trans only',
                     ifelse(reg == 'cis+trans',
                     ifelse((prop.h<prop.p & prop.h>=.5)|(prop.h>prop.p & prop.h<=.5), 'cis+trans',
                     ifelse((prop.h<prop.p & prop.h<.5 & prop.p>.5), "trans Mo17",
                     ifelse((prop.h>prop.p & prop.h>.5 & prop.p<.5), "trans B73",
                     ifelse((prop.h<prop.p & prop.p<.5), "cis Mo17",
                     ifelse((prop.h>prop.p & prop.p>.5), "cis B73", "unc"))))),
                     'unc')))) %>%
        mutate(reg = factor(reg, levels = reg_levels)) %>%
        dplyr::select(gid,B73,Mo17,ref,alt,cft,prop.p,prop.h,reg)
        t_ase
    #}}}
}
