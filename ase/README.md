# Pipeline to Analyze Allele-Specific Expression

## Requirements
* software
  * python3
  * an aligner: [star](https://github.com/alexdobin/STAR) or [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
  * [bedtools](https://bedtools.readthedocs.io/en/latest/)
  * custom python script: [ase.py](https://github.com/orionzhou/maize/blob/master/projects/ase.py)
* input files
  * [B73 AGPv4 genome fasta](https://www.maizegdb.org/assembly)
  * [AGPv4 gene interval file]()
  * [Mo17 SNP VCF](https://conservancy.umn.edu/handle/11299/198135?show=full)

## Align RNA-Seq reads to the reference
* (You can skip this step if you have already mapped your reads);
* The easiest way is to align everything to a single reference (B73 AGPv4) using [STAR](https://github.com/alexdobin/STAR);
* But to avoid mapping bias (e.g., some Mo17 reads contain mismatches (SNPs) to the B73 reference, so in general Mo17 samples will have a lower mapping rate than B73 samples, this is not a biological but technical bias), it is recommended to use [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) to build a variant-aware genome index using [Mo17 SNPs](https://conservancy.umn.edu/handle/11299/198135?show=full);
* [optional] Proper post-processing of obtained BAM alignments (PCR duplicate removal, etc.)

## Counting allele-specific reads
* suppose we have the RNA-Seq reads of a F1 hybrid mapped to B73 and obtained a BAM file called `01.bam`
* we also have the Mo17 SNPs downloaded and renamed to `Mo17.vcf.gz`
* the gene interval file (in BED format) is also downloaded and renamed to `genes.bed`
* create a folder and move all three files in:

      mkdir ase_test
      cd ase_tes
      mv ???/01.bam ????/Mo17.vcf.gz ?????/genes.bed ./

* convert BAM file to interval BED file and sort (this may take a while depending on the size of the BAM):

      ase.py bam2bed 01.bam 02.bed
      sort -k1,1 -k2,2n 02.bed -o 03.sorted.bed

* extract reads overlapping known Mo17 SNPs and sort:

      bedtools intersect -wa -wb -a 03.sorted.bed -b Mo17.bed > 04.ovlp.bed
      sort -k4,4 -k1,1 -k2,2n 04.ovlp.bed > 05.sorted.bed

* prepare extracted reads for B73/Mo17 allele separation and sort:

    ase.py bed_prep 05.sorted.bed 06.tsv 06.bed
    sort -k1,1 -k2,2n 06.bed -o 06.sorted.bed

* extract reads overlapping gene intervals:

    bedtools intersect -wa -wb -a genes.bed -b 06.sorted.bed > 07.bed

* summarise allele-specific reads at the gene level:

    ase.py bed_summarise 06.tsv 07.bed 10.tsv

* [optional] remove all intermediate files:

    rm 0[1-7].*

The final output `10.tsv` contains gene-level ASE reads counts.
