#!/bin/bash -l
#$ -N test
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -q all.q
#$ -l h_vmem=20G
#$ -l s_rt=10:00:00
#$ -pe smp 4
#$ -m a
#$ -M myemail@caas.cn
#$ -j y
#$ -o job.out

blastall -p blastn -i test.fa -d wheat.db -o abc.out
