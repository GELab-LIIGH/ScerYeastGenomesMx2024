python scripts/make_bed.py -g data/ref/SACE_S288C_v1_allChr.gff -i list_ortho.txt -o data/ortholog_regions.bed
snakemake -s Snakefile  -p -j 30 --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' all
