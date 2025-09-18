mkdir -p data/bams/ 
mkdir -p data/stats/
mkdir -p tmp
mkdir -p st

bash stage1.sh
cd ./backbone/
bash run.sh
Rscripts graph.R data/introgressions/ phylogeny.csv trueintros
