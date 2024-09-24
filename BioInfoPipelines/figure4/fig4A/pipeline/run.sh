mkdir -p data/
mkdir -p data/phyl
mkdir -p data/vcfs
mkdir -p st/
snakemake -s MakeBlocks.smk -p -j 30 --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' data/summary_statistics_blocks.csv
snakemake -s coords.smk -p -j 10 --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' coords_files
snakemake -s Snakefile -p -j 30 --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' data/done.txt
