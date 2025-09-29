# de ArtSACE/Admixture_Fig1/Journal


qsub Remove_LD_from_VCF_SACE467.sge  /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/Matrix_SNPs_SACE_from_CONC_gt_onlySNPs_filtered_missing_10_biallelic 50 5 .5 
qlogin
cd /mnt/Timina/lmorales/Public/ymez/bin/scripts/10_admixture
module load python37/3.7.0 
python3.7 SendAdmixtureSGEs.py --dirbed /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/Admixture --bed Matrix_LDPrune_50_5_.5.bed --out /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/Admixture/out --minS 1 --maxS 5 --minK 1 --maxK 35 --setname SACE467

mv javelar_admixture_K1-35_S1-5.sh javelar_admixture_K1-35_S1-5_SACE467.sh 
bash javelar_admixture_K1-35_S1-5_SACE467.sh

cd /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/Admixture/out
grep -h CV /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE467/Admixture/out/log*.out|awk -F "=|:" '{print $2"\t"$3}' | sed 's/)//g'| sort -n #> Values_CV_SACE467.txt


/usr/bin/pong -m Pong_mfile_SACE467_K1-K10_SBest.txt -i ind2pop_SACE467_AllStrains.csv -p 5000 -n order_strains.csv  --color_list /mnt/Timina/lmorales/Public/ymez/tmp/06_genotyping/trees/SACE421_Admixture1/out/colores_pong

