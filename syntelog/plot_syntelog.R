#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))

ps <- ArgumentParser(description = 'Make synteny plot using provided fasta sequences and mark features')
ps$add_argument("f_fa", help="multi-fasta sequences")
ps$add_argument("fp", help="output (PDF) file")

ps$add_argument("--gid", default='test gene', help="gene ID used to extract structure info [default: %(default)s]")
ps$add_argument("--bed", default=NULL, help="BED file containing feature positions to highlight")
ps$add_argument("--gene", help="R data file containing gene structure information for provided orthologs [default: %(default)s]",
    default="/home/springer/zhoux379/projects/genome/data2/syntelog/xref.maize.v4.rds")
ps$add_argument("--width", type='integer', help="figure width [default: %(default)s]", default=7)
ps$add_argument("--height", type='integer',  help="figure height [default: %(default)s]", default=6)
args <- ps$parse_args()

f_fa = args$f_fa
f_bed = args$f_bed
fp = args$fp
gid = args$gid; f_gene = args$gene
wd = args$width; ht = args$height

#{{{ read & functions
require(Biostrings)
require(tidyverse)
require(glue)
require(ggtree)
require(tidytree)
require(ggpubr)
otheme <- function(margin = c(.5,.5,.5,.5),
                   strip.size = 8, strip.style = 'white', strip.compact = F,
                   legend.pos='', legend.dir='v', legend.box='v', legend.vjust=NA,
                   legend.title = F, legend.border = F,
                   legend.spacing.x = 0, legend.spacing.y = 0,
                   panel.border = T, panel.spacing = .02,
                   xticks = F, yticks = F, xtitle = F, ytitle = F,
                   xtext = F, ytext = F, xgrid = F, ygrid = F,
                   xsize = 8, ysize = 8, nostrip = T) {
    #{{{ custom theme
    o = theme_bw()
    #{{{ strip style
    strip.color = ifelse(strip.style=='dark', 'white', 'black')
    strip.fill = ifelse(strip.style=='dark', '#5D729D', 'white')
    strip.box.color = ifelse(strip.style=='dark', '#4A618C', '#4A618C')
    if(strip.style == 'dark') {
        o = o +
            theme(strip.background = element_blank(),
                  strip.text = element_textbox(size = strip.size,
                  color = strip.color, fill = strip.fill, box.color = strip.box.color,
                  halign=.5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
                  padding = margin(2, 0, 1, 0), margin = margin(.05, .05, .05, .05)))
    } else {
        o = o +
            theme(strip.background = element_blank(),
                  strip.text = element_text(size = strip.size))
        if(strip.compact)
            o = o + theme(strip.text=element_text(margin=margin(.05,.05,.05,.05)))
    }
    #}}}
    #{{{ legend
    o = o +
        theme(legend.background = element_blank(),
              legend.key.size = unit(.8, 'lines'),
              legend.text = element_text(size = 8),
              legend.spacing.x = unit(legend.spacing.x, "lines"),
              legend.spacing.y = unit(legend.spacing.y, "lines"),
              legend.margin = margin(.1,.1,.1,.1,'lines'))
    if(legend.title == F) {
        o = o + theme(legend.title = element_blank())
    } else {
        o = o + theme(legend.title = element_text(size=8))
    }
    if(legend.pos == 'none') {
        o = o + theme(legend.position = 'none')
    } else if(legend.pos == 'top.center.out') {
        if(is.na(legend.vjust)) legend.vjust=0
        o = o + theme(legend.position = c(.5,1), legend.justification = c(.5,legend.vjust))
        margin[1] = margin[1] + 1
    } else if(legend.pos == 'top.center') {
        o = o + theme(legend.position = c(.5,1), legend.justification = c(.5,1))
    } else if(legend.pos == 'top.left') {
        o = o + theme(legend.position = c(0,1), legend.justification = c(0,1))
    } else if(legend.pos == 'top.right') {
        o = o + theme(legend.position = c(1,1), legend.justification = c(1,1))
    } else if(legend.pos == 'bottom.right') {
        o = o + theme(legend.position = c(1,0), legend.justification = c(1,0))
    } else if(legend.pos == 'bottom.left') {
        o = o + theme(legend.position = c(0,0), legend.justification = c(0,0))
    } else if(legend.pos == 'left') {
        o = o + theme(legend.position = 'left')
    }
    if(legend.dir == 'h')
        o = o + theme(legend.direction='horizontal')
    if(legend.box == 'h')
        o = o + theme(legend.box='horizontal')
    if(legend.border == T)
        o = o + theme(legend.box.background = element_rect())
    #}}}
    #{{{ axis tick/text/title
    if(xtitle) {
        o = o + theme(axis.title.x = element_text(size = 9))
    } else {
        o = o + theme(axis.title.x = element_blank())
    }
    if(ytitle) {
        o = o + theme(axis.title.y = element_text(size = 9))
    } else {
        o = o + theme(axis.title.y = element_blank())
    }
    if(xtext) {
        o = o + theme(axis.text.x = element_text(size = xsize))
    } else {
        o = o + theme(axis.text.x = element_blank())
    }
    if(ytext) {
        o = o + theme(axis.text.y = element_text(size = ysize))
    } else {
        o = o + theme(axis.text.y = element_blank())
    }
    if(!xticks)
        o = o + theme(axis.ticks.x = element_blank(), axis.ticks.length.x=unit(0,'pt'))
    if(!yticks)
        o = o + theme(axis.ticks.y = element_blank(), axis.ticks.length.y=unit(0,'pt'))
    #}}}
    #{{{ panel & grid
    o = o + theme(panel.grid.minor = element_blank()) +
        theme(panel.spacing = unit(panel.spacing, "lines"))
    if(!xgrid)
        o = o + theme(panel.grid.major.x = element_blank())
    if(!ygrid)
        o = o + theme(panel.grid.major.y = element_blank())
    if(!panel.border)
        o = o + theme(panel.border = element_blank())
    #}}}
    o + theme(plot.margin = unit(margin, "lines")) +
        theme(title=element_text(margin=margin(0,0,0,0,'cm')))
    #}}}
}
read_bed  <- function(f_bed) {
#{{{
    if (!is.null(f_bed) & file.exists(f_bed)) {
        bed = read_tsv(f_bed, col_names=c("sid",'beg','end'))
        if (str_detect(bed$sid[[1]], '%')) {
            bed = bed %>% separate('sid',c('mid','sid'), sep='%') %>% select(-mid)
        }
        bed %>% mutate(pos = (beg+end)/2) %>%
            select(gt=sid, beg, end, pos)
    } else {
        tibble()
    }
#}}}
}
msa_tree <- function(seqs) {
    #{{{
    pre = glue("tmp.msa.{sample(10000, 1)}")
    writeXStringSet(seqs, glue('{pre}.fa'))
    system(glue("muscle -in {pre}.fa -out {pre}.2.fa -tree2 {pre}.2.nwk"))
    msa = readDNAMultipleAlignment(filepath=glue("{pre}.2.fa"), format='fasta')
    tree = read.tree(glue("{pre}.2.nwk"))
    labs = with(subset(fortify(tree), isTip), label[order(y, decreasing=T)])
    nl = length(labs)
    ty = tibble(lab = labs, y = nl:1) %>%
        separate(lab, c('gt','extra'), sep='_', remove=F) %>%
        select(gt, lab, y)
    #
    system(glue("rm {pre}.*"))
    list(msa=msa, tree=tree, ty=ty)
    #}}}
}
make_syn  <- function(name1, name2, seqs) {
    #{{{ extract synteny blocks and mismatch positions from alignment
    require(Biostrings)
    #x = unmasked(msa)
    #s1 = RemoveGaps(x[name1])[[1]]
    #s2 = RemoveGaps(x[name2])[[1]]
    pw = pairwiseAlignment(seqs[[name1]], seqs[[name2]], type='local')
    toff = start(pattern(pw)); qoff = start(subject(pw))
    y1 = as.character(pattern(pw))
    y2 = as.character(subject(pw))
    h0 = as.character(compareStrings(y2,y1))
    h = str_replace_all(h0, "[ATCGN?]", '*')
    h2=unlist(strsplit(h, split = ""))
    j = tibble(v = rle(h2)$values, aSize = rle(h2)$lengths) %>%
        mutate(aEnd=cumsum(aSize), aBeg = aEnd-aSize) %>%
        mutate(tSize = ifelse(v=='+', 0, aSize)) %>%
        mutate(qSize = ifelse(v=='-', 0, aSize)) %>%
        mutate(tEnd=cumsum(tSize), tBeg=tEnd-tSize) %>%
        mutate(qEnd=cumsum(qSize), qBeg=qEnd-qSize) %>%
        select(v, aBeg,aEnd,aSize,tBeg,tEnd,tSize,qBeg,qEnd,qSize) %>%
        mutate(tBeg = toff+tBeg-1, tEnd=tBeg+tSize) %>%
        mutate(qBeg = qoff+qBeg-1, qEnd=qBeg+qSize)
    get_tpos <- function(v, aSize, tBeg, tEnd)
        ifelse(v=='+', list(rep(tBeg, aSize)), list((tBeg+1):tEnd))
    get_qpos <- function(v, aSize, qBeg, qEnd)
        ifelse(v=='-', list(rep(qBeg, aSize)), list((qBeg+1):qEnd))
    j2 = j %>%
        mutate(tpos = pmap(list(v, aSize, tBeg, tEnd), get_tpos)) %>%
        mutate(qpos = pmap(list(v, aSize, qBeg, qEnd), get_qpos))
    tposs = unlist(j2$tpos)
    qposs = unlist(j2$qpos)
    mm = str_locate_all(h0, "\\?")[[1]] %>% as_tibble() %>%
        mutate(qPos = qposs[start], tPos = tposs[start]) %>%
        select(aPos=start,tPos,qPos)
    list(aln = j, mm = mm)
    #}}}
}
plot_syn <- function(ty, seqs) {
    #{{{
    #{{{ prepare syn tible
    labs = ty %>% arrange(desc(y)) %>% pull(lab)
    to = tibble(tgt=labs[1:nrow(ty)-1], qry=labs[2:nrow(ty)]) %>%
        mutate(x = map2(tgt,qry, make_syn, seqs=seqs)) %>%
        mutate(syn = map(x, 'aln'), mm = map(x, 'mm')) %>%
        select(tgt, qry, syn, mm)
    #
    tp = to %>% select(tgt,qry,syn) %>% unnest(syn) %>%
        filter(v == '*') %>%
        inner_join(ty %>% select(lab,y), by=c('tgt'='lab')) %>% rename(y1=y) %>%
        inner_join(ty %>% select(lab,y), by=c('qry'='lab')) %>% rename(y2=y) %>%
        mutate(tBeg=tBeg + 1, qBeg = qBeg + 1) %>%
        mutate(i = 1:n())
    tp1 = tp %>% select(i, y=y1, tBeg, tEnd) %>%
        gather(type, pos, -i, -y) %>% mutate(y=y-.1)
    tp2 = tp %>% select(i, y=y2, qBeg, qEnd) %>%
        gather(type, pos, -i, -y) %>% mutate(y=y+.1)
    coordmap = c("tBeg"=1,'tEnd'=2,'qEnd'=3,'qBeg'=4)
    tp = tp1 %>% rbind(tp2) %>% mutate(i2 = coordmap[type]) %>% arrange(i, i2)
    #}}}
    #
    p_syn = ggplot(tp) +
        geom_polygon(aes(x=pos,y=y,group=i), fill='royalblue', alpha=.2,
                     size=0,color=NA) +
        coord_cartesian(xlim = c(0,4000)) +
        scale_x_continuous(breaks=c(0,2000,4000), labels=c('-2k','TSS','+2k'),
                           expand=expansion(mult=c(.01,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.01,.01))) +
        otheme(xtext=T,xtick=T,ytext=T,ytick=T, panel.border=F,
            margin=c(.2,.5,.2,0))
    p_syn
    #}}}
}
add_gene_track <- function(gid, tg, ty, p, ht.exon=.05, ht.cds=.1) {
    #{{{
    tg1 = tg %>% filter(gid==!!gid) %>% inner_join(ty, by='gt') %>%
        mutate(beg = beg+2000, end = end+2000)
    tgr = tg1 %>% group_by(gt, gid,y) %>%
        summarise(beg=min(beg), end=max(end)) %>% ungroup()
    tge = tg1 %>% filter(type=='exon')
    tgc = tg1 %>% filter(type=='CDS')
    col.intron='grey'; col.exon='grey'; col.cds='royalblue'
    arw = arrow(length=unit(.1,'cm'), angle=30, ends='last',type="open")
    p +
        geom_segment(data=tgr,aes(x=beg,xend=end,y=y,yend=y),col=col.intron,size=.5) +
        geom_rect(data=tge,aes(xmin=beg,xmax=end,ymin=y-ht.exon,ymax=y+ht.exon),fill=col.exon,color=NA,alpha=1) +
        geom_rect(data=tgc,aes(xmin=beg,xmax=end,ymin=y-ht.cds,ymax=y+ht.cds),fill=col.cds,color=NA,alpha=1) +
        geom_segment(data=tgr,aes(x=2000,xend=2100,y=y+ht.cds*1.1,yend=y+ht.cds*1.1),
                     color='black', size=.2, arrow=arw) +
        geom_rect(xmin=1999,xmax=2001,ymin=-Inf,ymax=Inf, fill='yellow', alpha=.2)
    #}}}
}
plot_title <- function(gid, gname=NA,gnote=NA, font.size=3) {
    #{{{ #plot title (motif + gene ID + alias + annotation)
    glab = ifelse(is.na(gname), gid, glue("{gid} ({gname})"))
    glab = ifelse(is.na(gnote), glue("{glab} (unknown)"), glue("{glab}: {gnote}"))
    #
    tit = glue("{glab}")
    p02 = ggplot() +
        annotate('text', x=0, y=0, label=tit, size=font.size, hjust=.5, vjust=.5) +
        otheme(panel.border=F, margin=c(.2,.2,0,.2))
    p02
    #}}}
}
combo_plot  <- function(gid, ty, seqs, bed, tree, tg) {
    #{{{
    p0 = plot_title(gid, font.size=3)
    # tree plot
    ty2 = ty %>% select(taxa=lab, gt)
    p_tree = ggtree(tree) %<+% ty2 +
        #geom_tiplab(aes(label=gt), hjust=0, align=T) +
        scale_x_continuous(expand=expansion(mult=c(.05,.005))) +
        scale_y_continuous(expand=expansion(mult=c(.01,.05))) +
        otheme(panel.border=F, margin=c(.2,.2,.2,0))
    #
    p1 = plot_syn(ty, seqs)
    mtf2 = bed %>% inner_join(ty, by='gt')
    #umr2 = umr %>% inner_join(ty, by='gt')
    pg = add_gene_track(gid,tg,ty,p1) +
        #geom_segment(data=umr2, aes(x=start,xend=end,y=y-.1,yend=y-.1), color='green',
                     #size=.5, alpha=1) +
        geom_point(data=mtf2,aes(x=pos, y=y-.1),col='red',size=1)
    #
    #p_lfc = plot_lfc(gid, ty, lfc)
    #
    #p_msa = plot_msa(aln, ty)
    #pb = ggarrange(p_tree, pg, p_msa, p_lfc, nrow=1, ncol=4, widths=c(1,5,.5,1))
    #pb = ggarrange(p_tree, pg, p_lfc, nrow=1, ncol=3, widths=c(1,5,1))
    pb = ggarrange(p_tree, pg, nrow=1, ncol=2, widths=c(1,5))
    ggarrange(p0, pb, nrow=2, heights=c(1,20))
    #}}}
}
#}}}
tg = readRDS(f_gene)

seqs = readDNAStringSet(f_fa)
x = msa_tree(seqs)
msa = x$msa; tree = x$tree; ty = x$ty
bed = read_bed(f_bed)
p = combo_plot(gid, ty, seqs, bed, tree, tg)
p %>% ggexport(filename=fp, width=wd, height=ht)


