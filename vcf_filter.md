## Filter a VCF file using various filters

This code uses [bcftools view](http://samtools.github.io/bcftools/bcftools.html#view) command to filter an input VCF file

Code:

    bcftools view -v snps -m2 -M2 -i 'F_PASS(GQ>=20 & GT!="mis") >= 0.3 && N_PASS(GT="AA") > 0 && MAF>=0.05' -Oz -o clean.vcf.gz All.final_snps.vcf.gz

This line of code applies 3 filters:
* `-v snps -m2 -M2`: only keep bi-allelic SNPs
* `F_PASS(GQ>=20) & GT!="mis") >= 0.3`: at least 30% of samples have non-missing genotype calls with genotype quality (`GQ`) higher than 20
* `MAF>=0.05`: minor allele frequency above 5%
