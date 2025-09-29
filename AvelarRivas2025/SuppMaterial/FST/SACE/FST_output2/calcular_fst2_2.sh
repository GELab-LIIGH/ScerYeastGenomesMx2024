#!/bin/bash
## Use current working directory (change working directory)
##Error file
#$ -e /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output/FST_SACE_$JOB_ID.err
## Out file
#$ -o /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output/FST_SACE_$JOB_ID.out
#$ -S /bin/bash
## Job's name
#$ -N GetFSTMatrix
#$ -l vf=8G
#$ -pe openmp 4
#$ -m e
## notification
#$ -M jabrahamavelar@gmail.com
## Modules
module load vcftools/0.1.14 
#!/bin/bash

VCF_FILE="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/Matrix_SNPs_SACE_from_SACE_gt_onlySNPs_filtered_missing_10_biallelic"
mkdir -p /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output2
# Define solo las poblaciones deseadas
POP_DIR="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/Grupos"
POP_FILES=(
  "$POP_DIR/MMMG.txt"
  "$POP_DIR/FG.txt"
  "$POP_DIR/Tamaulipas.txt"
  "$POP_DIR/SAM2_b.txt"
  "$POP_DIR/WB3.txt"
  "$POP_DIR/Tequila.txt"
)
SUMMARY_FILE="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output2/fst_summary.txt"
echo -e "Comparación\tMean_Fst\tWeighted_Fst" > "$SUMMARY_FILE"
# Iterar sobre todos los pares posibles de poblaciones (todas vs todas)

# ##!/bin/bash
#VCF_FILE="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/Matrix_SNPs_SACE_from_SACE_gt_onlySNPs_filtered_missing_10_biallelic"
#mkdir -p /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output
#POP_FILES=($(ls /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/Grupos/*.txt))
#SUMMARY_FILE="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output/fst_summary.txt"
echo -e "Comparación\tMean_Fst\tWeighted_Fst" > "$SUMMARY_FILE"

# Iterar sobre todos los pares de poblaciones
for ((i=0; i<${#POP_FILES[@]}; i++)); do
    for ((j=i+1; j<${#POP_FILES[@]}; j++)); do
        POP1=${POP_FILES[i]}
        POP2=${POP_FILES[j]}
        
        # Obtener nombres de las poblaciones sin la ruta y la extensión
        NAME1=$(basename "$POP1" .txt)
        NAME2=$(basename "$POP2" .txt)

        # Definir nombre de salida
        OUTPUT="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/FST_output/fst_${NAME1}_vs_${NAME2}"

        echo "Calculando Fst para $NAME1 vs $NAME2..."

        # Ejecutar vcftools para calcular Fst
        vcftools --gzvcf "$VCF_FILE" --weir-fst-pop "$POP1" --weir-fst-pop "$POP2" --out "$OUTPUT" > "${OUTPUT}.log"

        # Extraer valores de Fst del log de vcftools
        MEAN_FST=$(grep "Weir and Cockerham mean Fst estimate" "${OUTPUT}.log" | awk '{print $6}')
        WEIGHTED_FST=$(grep "Weir and Cockerham weighted Fst estimate" "${OUTPUT}.log" | awk '{print $7}')

        # Guardar en el archivo de resumen
        echo -e "${NAME1}_vs_${NAME2}\t${MEAN_FST}\t${WEIGHTED_FST}" >> "$SUMMARY_FILE"
    done
done

echo "Análisis completado. Resultados en /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE489/SACE/ $SUMMARY_FILE"
