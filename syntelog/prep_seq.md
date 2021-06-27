## Prepare fasta sequences for synteny visualization

### Configuration files:
- Genome configure file (`genomes.yml`) containing paths to genome fasta ang gene GFF files, etc
  - available here: `/datalus/weiyu/projects/genome/nf/genomes.yml`
- R data file containing gene structure information for provided orthologs
  - pre-made file available for:
    - maize B73v4-anchored ortholog structure in 30 maize genotypes: `/datalus/weiyu/projects/genome/data2/syntelog/xref.maize.v4.rds`
    - maize B73v5-anchored ortholog structure in 30 maize genotypes: `/datalus/weiyu/projects/genome/data2/syntelog/xref.maize.v5.rds`
    - wheat CS-D-anchored ortholog structure in 30 maize genotypes: `/datalus/weiyu/projects/genome/data2/syntelog/xref.wheat.csd.rds`

Run the script:

    ./prep_seq.R --opt genomic+2k Zm00001d003460 test1.fa

For detailed help:

```bash
(wgc) pzhou@sugon1:syntelog $ ./prep_seq.R -h
usage: ./prep_seq.R [-h] [--opt OPT] [--gene GENE] [--gcfg GCFG] gid out

Prepare fasta sequences of syntelogs

positional arguments:
  gid          gene ID used to extract structure info
  out          output multi-fasta file

optional arguments:
  -h, --help   show this help message and exit
  --opt OPT    range options [default: cds]
  --gene GENE  R data file containing gene structure information for provided
               orthologs [default: /datalus/weiyu/projects/genome/data2/syntel
               og/xref.maize.v4.rds]
  --gcfg GCFG  genome configure file (yml) [default:
               /datalus/weiyu/projects/genome/nf/genomes.yml]
```
