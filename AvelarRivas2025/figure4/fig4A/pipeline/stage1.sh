######################################
##                                  ##
##Automatization to obtain figure 4A##
##                                  ##
######################################

######################################
#Constants#
sge_dir="st"
st_dir=$sge_dir
bam_dir="data/bams/"
stats_dir="data/stats/"
tmp_dir="tmp/"
script_dir="scripts/"
gff_folder="data/"
data_folder="data/"
intro_dir=$data_folder"introgressions/"
ortho_file=$data_folder"41588_2017_BFng3847_MOESM339_ESM.tsv"
sace_gff_file=$gff_folder"SACE_S288C_v1_allChr.gff"
sapa_gff_file=$gff_folder"SAPA_YPS138_v1.gff"
sace_txt_file=$data_folder"SACE_S288C_v1_allChr.txt"
sapa_txt_file=$data_folder"SAPA_YPS138_v1.txt"
launch_cv_stats_sges=$sge_dir"coverage_plot.sh"
cv_stats_sge_suffix="_coverage_stats.sge"
launch_ov_genes_sges=$sge_dir"overlapping_genes.sh"
ov_genes_sge_suffix="_overlapping_genes.sge"
launch_introgression_sges=$sge_dir"introgressions.sh"
introgression_sge_suffix="_introgression.sge"
######################################

##
module load python37/3.7.0 htslib/1.9 bedtools/2.24
##

echo "Current wordking directory:"
echo "$PWD"


echo " "
echo "========================================"
echo "    Writing coverage statistics sges    "
echo "========================================"
echo " "

	echo "#Launches al sges" > $launch_cv_stats_sges
for sample in $(ls data/bams/*_CONC.rmdup.addgp.bam|cut -d"_" -f1|cut -d"/" -f3)
do
	sge_file=$sge_dir$sample$cv_stats_sge_suffix
	bam_file=$bam_dir$sample"_CONC.rmdup.addgp.bam"
	depth_file=$tmp_dir$sample"_CONC_Q20.depth"
	echo "Writing "$sge_file" ..."
	echo "#!/bin/bash" > $sge_file
	echo "#$ -cwd" >> $sge_file
	echo ". /etc/profile.d/modules.sh" >> $sge_file
	echo "#$ -e " $st_dir$sample"_cv_stats.err" >> $sge_file
	echo "#$ -o " $st_dir$sample"_cv_stats.out" >> $sge_file
	echo "#$ -S /bin/bash" >> $sge_file
	echo "#$ -N " $sample"cv_stats" >> $sge_file
	echo "#$ -l vf=4G" >> $sge_file
	echo "#$ -pe openmp 10" >> $sge_file
	echo "source /etc/bashrc" >> $sge_file
	echo "module load htslib/1.2.1 gcc/5.1.0 samtools/1.9 r/3.6.1" >> $sge_file
	echo "samtools depth -aa -Q 20 "$bam_file" > "$depth_file >> $sge_file
	echo "R CMD BATCH --no-save --no-restore \"--args FILE='"$depth_file"' SIZE=1000 STEP=1000 SAMPLE='"$sample"' STATS='"$stats_dir"'\" "$script_dir"make_depth_plot.R "$tmp_dir$sample"_make_depth_plot.Rout" >> $sge_file
	echo "rm "$depth_file >> $sge_file
	echo "qsub " $sge_file >> $launch_cv_stats_sges
	echo "Done!"
done

echo " "
echo "========================================"
echo "Obtaining coverage statistics"
echo "========================================"
echo " "

bash $launch_cv_stats_sges
while [[ $(qstat |wc -l) > 0 ]];
do
	echo "Jobs remaining: "$((x=$(qstat | wc -l), y=2, x-y))
	sleep 1m
done

echo "Done!"

echo " "
echo "========================================"
echo "     Converting gffs to txt and bed     "
echo "========================================"
echo " "

python3 ${script_dir}gff2txt.py -g $sapa_gff_file -f $data_folder
python3 ${script_dir}gff2txt.py -g $sace_gff_file -f $data_folder


echo "Done!"

echo " "
echo "========================================"
echo "     Writing overlapping genes sges     "
echo "========================================"
echo " "


echo "#Launches al sges" > $launch_ov_genes_sges
for sample in $(ls data/bams/*_CONC.rmdup.addgp.bam|cut -d"_" -f1|cut -d"/" -f3)
do
	sge_file=$sge_dir$sample$ov_genes_sge_suffix
	bam_file=$bam_dir$sample"_CONC.rmdup.addgp.bam"
	windows_file=$stats_dir$sample"_windows.csv"
	echo "Writing "$sge_file" ..."
	echo "#!/bin/bash" > $sge_file
	echo "#$ -cwd" >> $sge_file
	echo ". /etc/profile.d/modules.sh" >> $sge_file
	echo "#$ -e " $st_dir$sample"_overlapping_genes.err" >> $sge_file
	echo "#$ -o " $st_dir$sample"_overlapping_genes.out" >> $sge_file
	echo "#$ -S /bin/bash" >> $sge_file
	echo "#$ -N " $sample"_ov_genes" >> $sge_file
	echo "#$ -l vf=4G" >> $sge_file
	echo "#$ -pe openmp 2" >> $sge_file
	echo "source /etc/bashrc" >> $sge_file
	echo "module load htslib/1.9 gcc/5.1.0 bedtools/2.26 python37/3.7.0 samtools/1.9" >> $sge_file
	echo "python3 "$script_dir"overlapping_genesAbove.py -t "$sapa_txt_file" -f "$stats_dir" -w "$windows_file" -b "$bam_file >> $sge_file
	echo "python3 "$script_dir"overlapping_genesAbove.py -t "$sace_txt_file" -f "$stats_dir" -w "$windows_file" -b "$bam_file >> $sge_file
	echo "qsub " $sge_file >> $launch_ov_genes_sges
	echo "Done!"
done


echo " "
echo "========================================"
echo "      Obtaining overlapping genes       "
echo "========================================"
echo " "

bash $launch_ov_genes_sges
while [[ $(qstat |wc -l) > 0 ]];
do
	echo "Jobs remaining: "$((x=$(qstat | wc -l), y=2, x-y))
	sleep 1m
done

echo "Done!"

echo " "
echo "========================================"
echo "       Writing introgression sges       "
echo "========================================"
echo " "


echo "#Launches al sges" > $launch_introgression_sges
for sample in $(ls data/bams/*_CONC.rmdup.addgp.bam|cut -d"_" -f1|cut -d"/" -f3)
do
	sge_file=$sge_dir$sample$introgression_sge_suffix
	bam_file=$bam_dir$sample"_CONC.rmdup.addgp.bam"
	sace_file=$stats_dir$sample"_SACE_genes_50aboveQm.csv"
	sapa_file=$stats_dir$sample"_SAPA_genes_50aboveQm.csv"
	perRef_file=$stats_dir$sample"_CONC_Q20_perReference.csv"
	echo "Writing "$sge_file" ..."
	echo "#!/bin/bash" > $sge_file
	echo "#$ -cwd" >> $sge_file
	echo ". /etc/profile.d/modules.sh" >> $sge_file
	echo "#$ -e " $st_dir$sample"_introgression.err" >> $sge_file
	echo "#$ -o " $st_dir$sample"_introgression.out" >> $sge_file
	echo "#$ -S /bin/bash" >> $sge_file
	echo "#$ -N " $sample"_ov_genes" >> $sge_file
	echo "#$ -l vf=4G" >> $sge_file
	echo "#$ -pe openmp 2" >> $sge_file
	echo "source /etc/bashrc" >> $sge_file
	echo "module load htslib/1.9 gcc/5.1.0 bedtools/2.26 python37/3.7.0 samtools/1.9" >> $sge_file
	echo "python3 "$script_dir"ortholog_coverageAll.py -s "$sample" -f "$data_folder" -a "$sapa_file" -e "$sace_file" -b "$bam_file" -O "$intro_dir" -o "$ortho_file >> $sge_file
	echo "qsub " $sge_file >> $launch_introgression_sges
	echo "Done!"
done

echo " "
echo "========================================"
echo "      Obtaining introgressed genes      "
echo "========================================"
echo " "

bash $launch_introgression_sges
while [[ $(qstat |wc -l) > 0 ]];
do
	echo "Jobs remaining: "$((x=$(qstat | wc -l), y=2, x-y))
	sleep 1m
done

echo "Done!"
