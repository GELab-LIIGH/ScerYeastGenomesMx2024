mkdir -p data/bams/ 
mkdir -p data/stats/
mkdir -p tmp
mkdir -p st

cd data/
tar -zxf CONC_Mez1_v1_allChr.fasta.tar.gz
tar -zxf CONC_Mez1_v1_allChr.fasta.sa.tar.gz

bash stage1.sh
cd trueIntrogressions/
bash run.sh
Rscripts graph.R data/introgressions/ phylogeny.csv trueintrosFile
