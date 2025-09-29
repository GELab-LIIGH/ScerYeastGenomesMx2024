module load vcftools/0.1.14
matriz="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/Matrix_SNPs_SACE_from_CONC_gt_onlySNPs_filtered_missing_10_biallelic" #../../Admixture/Matrix_SNPs_SACE_from_SACE_gt_onlySNPs_filtered_missing_10_biallelic  
matrizAll="/mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/Matrix_SNPs_SACE_from_CONC_gt.vcf"  #../../Admixture/Matrix_SNPs_SACE_from_SACE_gt.vcf 

for popfile in *.txt; # puede ser también csv; 
do 
echo $popfile;     
population="${popfile%.txt}"; 
echo $population;  

#vcftools --vcf "$matriz"  --keep "$popfile" --site-pi --out "${population}_results";
#vcftools --vcf "$matriz" --keep "$popfile"  --TajimaD 10000 --out "${population}_TDw_results";
#vcftools --vcf "$matriz"  --keep "$popfile" --window-pi 10000 --out "${population}_Window_results";
#vcftools --vcf "$matrizAll"  --keep "$popfile" --site-pi --out "${population}_results_AllMatrix";


vcftools --vcf "$matrizAll" --keep "$popfile" --window-pi 10000 --out "${population}_Window_results_AllMatrix";
vcftools --vcf "$matrizAll" --keep "$popfile"  --TajimaD 10000 --out "${population}_TDw_results_AllMatrix";

done