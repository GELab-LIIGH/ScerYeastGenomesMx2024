library(data.table)
library(ggplot2)

argsVal <- commandArgs(trailingOnly = T)
introDir <- argsVal[1]
samples <- argsVal[2]
out <- argsVal[3]

introFiles <- list.files(introDir,pattern = "_introSAPA_All_noOrth.csv")

strainV <- c()
noIntrosV <- c()
heterozygosityV <- c()
for(inF in introFiles){
  f <- fread(paste(introDir,inF,sep="/"))
  str <- strsplit(inF,split="_")[[1]][1]
  no <- dim(f)[1]
  strainV <- c(strainV,str)
  noIntrosV <- c(noIntrosV,no)
  heterozygosityV <- c(heterozygosityV,(dim(f[Heterozygosity == "Het"])[1])/no)
}
data <- data.frame(ID = strainV,noIntros = noIntrosV,heterozygosity = heterozygosityV)

smpData <- fread(samples)
data <- merge(data,smpData)
setDT(data)
data[PhyloGroup_SACE469 !="MMMG" & PhyloGroup_SACE469!="Tamaulipas" & PhyloGroup_SACE469!="Tequila" & PhyloGroup_SACE469!= "SAM2" & PhyloGroup_SACE469!="FG" & PhyloGroup_SACE469!="Alpechin"]$PhyloGroup_SACE469 <- "Others"

data$PhyloGroup_SACE469 <- factor(data$PhyloGroup_SACE469, levels <- c("Others","Alpechin","FG","SAM2","Tequila","MMMG","Tamaulipas"))

p <- ggplot(data,aes(noIntros,PhyloGroup_SACE469)) +
  geom_violin() + theme_classic() +
  geom_jitter(height=0.2,aes(color=heterozygosity),size=1.8,alpha=0.7,shape=16) +
  scale_color_gradient(low= "yellow",high= "red3")

ggsave(out,plot=p,width = 9, height = 6)
