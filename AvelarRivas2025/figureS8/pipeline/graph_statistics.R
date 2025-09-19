library(data.table)
library(ggplot2)
####INS####
file = "../figure4/fig4A/pipeline/trueIntrogressions/data/summary_statistics_blocks.csv"
dos17 = "../figure4/fig4A/pipeline/backbone/SampleSheet_217_Figure2.csv"
sace = "../figure4/fig4A/pipeline/phylogeny.csv"
het = "../figure4/fig4A/pipeline/trueIntrogressions/data/introgressions/hetero.csv"
##########
###OUTS###
blocksVSintro_pdf.3 = "./avgsize_intros_individual.pdf"
##########

sace_df = fread(sace)[,c("ID","PhyloGroup_SACE469")]
het_df = fread(het)
file_df = merge(sace_df,het_df,by="ID")

wb3 <- c("BR021c1","BR010c1","BR011c1")
file_df[file_df$ID %in% wb3,]$PhyloGroup_SACE469 <- "WB3"
file_df$PhyloGroup_SACE469 = factor(file_df$PhyloGroup_SACE469,levels=c("Tamaulipas","MMMG","Tequila","FG","Alpechin","SAM2","WB3","NortAmOak","Wine14","MixOrigin8"))
file_df = file_df[order(file_df$PhyloGroup_SACE469)]

###Distribution
file_df[file_df$PhyloGroup_SACE469 == "Tamaulipas",]$PhyloGroup_SACE469 <- "Mexican Agave 1"
file_df[file_df$PhyloGroup_SACE469 != "Mexican Agave 1" & file_df$PhyloGroup_SACE469 != "MMMG"& file_df$PhyloGroup_SACE469 != "Tequila" & file_df$PhyloGroup_SACE469 != "FG" & file_df$PhyloGroup_SACE469 != "Alpechin" & file_df$PhyloGroup_SACE469 != "SAM2" & file_df$PhyloGroup_SACE469 != "WB3",]$PhyloGroup_SACE469 <- "Others"
file_df[file_df$PhyloGroup_SACE469 == "MMMG",]$PhyloGroup_SACE469 <- "Mexican Agave 2"
file_df[is.na(file_df$PhyloGroup_SACE469)]$PhyloGroup_SACE469 <- "Others"
file_df$PhyloGroup_SACE469 <- factor(file_df$PhyloGroup_SACE469,levels = c("Others","Alpechin","FG","SAM2","WB3","Tequila","Mexican Agave 1","Mexican Agave 2"))

p7.3 <- ggplot(file_df,aes(x=meanGperB,y=noBlocks)) +
  geom_point(aes(color=heterozygosity2),size=1,alpha=0.6,shape=16) +
  theme_bw() + scale_color_gradient(low= "yellow",high= "red2") +
  theme(axis.text = element_text(size=9), axis.title =element_text(size=12)) +
  ylab("Number of introgression blocks") + xlab("Average number of genes per block") +
  facet_wrap(PhyloGroup_SACE469 ~ .,ncol=4)
ggsave(blocksVSintro_pdf.3,plot=p7.3,height=11,width=20,units="cm",dpi=300)
