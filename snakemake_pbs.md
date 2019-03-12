# Run snakemake through PBS job server

1. add these job properties to the `params` section of each rule:

      rule:
          input:  ...
          output: ...
          params:
              q: small
              nodes: 1
              ppn: 1
              runtime: 2:00:00
              mem: 8gb
              M: XXX@umn.edu
              m: a
              r: n
              N: jobname
              o: jobname.out
              e: jobname.err
          threads: 4
          shell: ...

notice that the number of processors (CPUs) is set by the `threads` parameter.

2. run snakemake like this:

      snakemake --cluster "qsub -q {params.q} -l nodes={params.nodes}:ppn={threads},mem={params.mem},walltime={params.runtime} -N {params.N} -M {params.M} -m {params.m} -r {params.r} -o {params.o} -e {params.e}"

or add this line to your `~/.bashrc`:

      alias snq='snakemake --cluster "qsub -q {params.q} -l nodes={params.nodes}:ppn={threads},mem={params.mem},walltime={params.runtime} -N {params.N} -M {params.M} -m {params.m} -r {params.r} -o {params.o} -e {params.e}"'

and then run:

      snq

PS: don't forget to do a dry run `snakemake -n` before starting the cluster job!
