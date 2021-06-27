#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))

ps <- ArgumentParser(description = 'Prepare fasta sequences of syntelogs')
ps$add_argument("gid", help="gene ID used to extract structure info")
ps$add_argument("out", help="output multi-fasta file")

ps$add_argument("--opt", help="range options [default: %(default)s]",
    default="cds")
ps$add_argument("--gene", help="R data file containing gene structure information for provided orthologs [default: %(default)s]",
    default="/datalus/weiyu/projects/genome/data2/syntelog/xref.maize.v4.rds")
ps$add_argument("--gcfg", help="genome configure file (yml) [default: %(default)s]", 
    default="/datalus/weiyu/projects/genome/nf/genomes.yml")
args <- ps$parse_args()

gid = args$gene

require(yaml)
require(tidyverse)
require(glue)

tg = readRDS(args$gene)
cfg = read_yaml(args$gcfg)
opt = args$opt

prep_coord <- function(ti, opt='cds') {
    #{{{
    srds = unique(ti$srd)
    stopifnot(length(srds) == 1)
    srd = srds[1]
    if (opt == 'cds') {
        tb = ti %>% filter(etype=='CDS') %>%
            select(chrom, start, end, srd)
    } else if (opt == 'genomic') {
        tb = ti %>% group_by(chrom, srd) %>%
            summarise(start = min(start), end = max(end)) %>%
            select(chrom, start, end, srd)
    } else if (opt == 'genomic+2k') {
        tb = ti %>% group_by(chrom, srd) %>%
            summarise(start = min(start), end = max(end)) %>%
            mutate(start = start-2000, end = end+2000) %>%
            select(chrom, start, end, srd)
    } else if (opt == 'genomic+1k') {
        tb = ti %>% group_by(chrom, srd) %>%
            summarise(start = min(start), end = max(end)) %>%
            mutate(start = start-1000, end = end+1000) %>%
            select(chrom, start, end, srd)
    } else if (opt == 'tss+2k') {
        tb = ti %>% group_by(chrom, srd) %>%
            summarise(tss = ifelse(srd=='-', max(end), min(start))) %>% ungroup() %>%
            mutate(tss = tss - 1) %>%
            mutate(start = tss-2000, end = tss+2000) %>% mutate(start = start + 1) %>%
            select(chrom, start, end, srd)
    } else {
        stop(glue("unknown seq option: {opt}\n"))
    }
    if (srd == '-') {
        tb %>% arrange(chrom, desc(start))
    } else {
        tb %>% arrange(chrom, start)
    }
    #}}}
}
extract_fa  <- function(tb, cfg, gt='Zmays_B73', opt='cds', id='test') {
    #{{{
    ti = prep_coord(tb, opt) %>%
        mutate(start=start-1, score='.', fid=1:n()) %>%
        select(chrom, start, end, fid, score, srd)
    pre = glue("tmp.ps.{sample(10000, 1)}")
    fb = glue("{pre}.bed")
    write_tsv(ti, fb, col_names=F)
    #
    db = cfg$genomes[[gt]]$fasta
    fs = glue("{pre}.tsv")
    cmd = glue("bedtools getfasta -tab -nameOnly -fi {db} -bed {fb} -fo {fs}")
    system(cmd)
    to = read_tsv(fs, col_names=c('fid','seq')) %>% arrange(fid)
    system(glue("rm -rf {pre}.*"))
    str_c(to$seq, collapse='', sep='')
    #}}}
}

to = tg %>% filter(gid==args$gid) %>%
    mutate(seq = pmap_chr(list(tb=gene, gt=gt, id=gt), extract_fa, cfg=!!cfg, opt=!!opt))
to1 = to %>% select(gt, seq) %>% mutate(gt = glue(">{gt}"))

write_delim(to1, args$out, col_names=F, delim="\n")
