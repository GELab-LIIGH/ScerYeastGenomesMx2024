mkdir -p data/
mkdir -p data/phyl
mkdir -p data/vcfs
mkdir -p st/
rm -f copy_vcfs.sh
awk -F"," '{gsub(/"/, "", $1);system("echo ln -s /mnt/Timina/lmorales/Public/ymez/tmp/05_vcalling/"$1"_CONC.g.vcf.gz data/ >> copy_vcfs.sh")}' ./SACE_all487.csv
rm -f copy_2sapa_vcfs.sh
awk -F"," '{gsub(/"/, "", $1);system("echo ln -s /mnt/Timina/lmorales/Public/ymez/tmp/05_vcalling/"$1"_SAPA.g.vcf data/ >> copy_2sapa_vcfs.sh")}' ./forced2sapa.csv
rm -f copy_sapa2conc_vcfs.sh
awk -F"," '{gsub(/"/, "", $1);system("echo ln -s /mnt/Timina/lmorales/Public/ymez/tmp/05_vcalling/"$1"_CONC.g.vcf.gz data/ >> copy_sapa2conc_vcfs.sh")}' ./SAPA_analysis.csv
bash copy_vcfs.sh
bash copy_2sapa_vcfs.sh
bash copy_sapa2conc_vcfs.sh
snakemake -s MakeBlocks.smk -p -j 30 --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' data/summary_statistics_blocks.csv
snakemake -s coords.smk -p -j 10 --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' coords_files
snakemake --keep-going -s backbone.smk -p -j 30 --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' all
bash scripts/selectModels.sh
snakemake --keep-going -s Snakefile -p -j 100 --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' all
bash scripts/toRelaunch.sh
Rscript scripts/trueIntrogressions.R
