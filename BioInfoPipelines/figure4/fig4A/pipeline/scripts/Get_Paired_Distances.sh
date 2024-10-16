module load bcftools/1.9

cd  data/vcfs/
genes=$(ls -d */)
genes_array=($genes)
${#genes_array[@]} # tamaño del arreglo

for ((i=0; i<${#genes_array[@]}; i++)); do
  echo "Processing directory: ${genes_array[$i]}"
  dir=${genes_array[$i]%/}
  gen="${dir%/}"
  echo $gen
  echo $i
  bcftools query -l $gen/Matrix_"$gen"_combinegvcfs_genotype_onlySNPs_filterlow_filteredSNPs_missing40_biallelic.vcf  > /cepas_$gen.txt

while IFS= read -r strain1; do
    current_pos=$(grep -n "^$strain1$" /cepas_$gen.txt | cut -d: -f1)  # Almacenar la posición actual en el archivo de cepas
    echo $current_pos
    sed "1,${current_pos}d" /cepas_$gen.txt  | while IFS= read -r strain2; do     # Leer las cepas restantes en el archivo de cepas
      diferencias=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -s "$strain1","$strain2" $gen/Matrix_"$gen"_combinegvcfs_genotype_onlySNPs_filterlow_filteredSNPs_missing40_biallelic.vcf | awk 'BEGIN { FS = "\t" } { split($5, g1, ":"); split($6, g2, ":"); diff = 0; for (i=1;i<=length(g1); i++) { if (g1[i] != g2[i]){ diff++; } } print diff }'  | grep 1 | wc -l)       # Calcular las diferencias entre las dos cepas
      echo "$strain1,$strain2,$diferencias" >> data/PairedDistances/paired_distances_sapa_$gen.csv       # Escribir el resultado en un archivo CSV
    done
  done < /cepas_$gen.txt
done
rm /cepas_$gen.txt
