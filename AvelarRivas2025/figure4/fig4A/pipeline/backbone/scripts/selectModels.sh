folder="data/phyl/backbone/"
prefix="Matrix_backbone_"
suffix="_filteredSNP_10missing_biallelic.phy.model.gz"
shared="data/shared_introgressions.txt"

echo "Gene,Model" > data/BestModels.csv
for gene in $(tail -n+2 $shared | cut -f1 | sort -u)
do
  modelsFile=$folder$prefix$gene$suffix
  if [[ -f $modelsFile ]]; then
    model=$( zcat $modelsFile | grep "best_model_BIC" | cut -d" " -f2 )
    echo $gene','$model >> data/BestModels.csv
  else
    echo $gene',MFP' >> data/BestModels.csv
  fi
done
