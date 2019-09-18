require(tidyverse)

fi = '/home/springer/nosha003/acr_te_acr_gene_pav_melt_na_cpm.txt'

# read data
ti = read_tsv(fi) %>%
    rename(PAV_status = value)

# ti2 has 5 columns: gid, Tissue, PAV_status, num_samples, FPKM (a list of numbers)
ti2 = ti %>% group_by(gid, Tissue, PAV_status) %>%
    summarise(num_samples = n(), FPKM = list(FPKM)) %>%
    ungroup()

# filter ti2 by requiring >=2 "Absent" samples and >=2 "Present" samples
ti3 = ti2 %>% select(gid, Tissue, PAV_status, num_samples) %>%
    spread(PAV_status, num_samples) %>%
    filter(Absent >= 2, Present >= 2)

# extract extract the gid | Tissue combinations from ti2
# ti4 now has four columns: gid, Tissue, Absent (FPKMs for the Absent class) and Present (FPKMs for the Present class)
ti4 = ti2 %>% inner_join(ti3[,c('gid','Tissue')], by = c('gid','Tissue')) %>%
    select(-num_samples) %>%
    spread(PAV_status, FPKM)

# do the ranksum test
# ti5 has 2 additional columns: test (a hypothesis test object) and pval (pvalue extracted from the test object)
ti5 = ti4 %>%
    mutate(test = map2(Absent, Present, wilcox.test)) %>%
    mutate(pval = map_dbl(test, 'p.value'))


# check here for explanation of the "map" function:
# https://purrr.tidyverse.org/reference/map.html


