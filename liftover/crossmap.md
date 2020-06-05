# create a temporary folder to work on

    $ mkdir test; cd test

# prepare your query BED file

    $ echo -e "M01\t1200000\t1201000" > M.bed
    $ cat M.bed
    M01     1200000 1201000

# create symbolic link to the chain file

    $ ln -sf /home/springer/zhoux379/projects/wgc/data/raw/Mo17_B73/10.B73_Mo17.chain BtoM.chain
    $ ln -sf /home/springer/zhoux379/projects/wgc/data/raw/Mo17_B73/10.Mo17_B73.chain MtoB.chain

# install [CrossMap](http://crossmap.sourceforge.net/#convert-bed-format-files) and run

    $ CrossMap.py bed MtoB.chain M.bed > B.bed
    @ 2020-06-05 16:20:14: Read the chain file:  MtoB.chain

    $ cat B.bed
    M01     1200000 1201000 ->      B01     794173  795173

Sometimes the result will be "Unmap" or "Split" depending on the chain mapping.


