dirw = '/home/springer/zhoux379/projects/misc/trevor_nhej'

diri = file.path(dirw, 'input_folder')
diro = file.path(dirw, 'output_folder')
pres = c(
    "dataset1.txt",
    "dataset2.txt",
    "dataset3.txt"
)

for (pre in pres) {
    # read input file and plot
    fi = sprintf("%s/%s", diri, pre)
    p = ggplot()
    # save plot to output file
    fo = sprintf("%s/%s.pdf", diri, pre)
    ggsave(p, file=fo, width=8, height=8)
}

