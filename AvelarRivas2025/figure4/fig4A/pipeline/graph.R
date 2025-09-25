library(data.table)
library(ggplot2)

argsVal <- commandArgs(trailingOnly = T)
introDir <- argsVal[1]
samples <- argsVal[2]
fList <- argsVal[3]
out <- paste0(introDir,"fig4A.pdf")
hetOut <- paste0(introDir,"hetero.csv")

fMatrix <- fread(fList)
fMatrix <- fMatrix[fMatrix$TrueIntrogressions,]
fMatrix <- fMatrix[fMatrix$TrueOrigin != "Sp_und_outgroup",]
fMatrix <- fMatrix[fMatrix$TrueOrigin != 'WINE-JS445c1SAPA',]

introFiles <- list.files(introDir,pattern = "_introSAPA_All_noOrth.csv")
strainV <- c()
noIntrosV <- c()
heterozygosityV <- c()
noBlocksV <- c()
for(inF in introFiles){
  f <- fread(paste(introDir,inF,sep="/"))
  str <- strsplit(inF,split="_")[[1]][1]
  tmpMat <- fMatrix[Strain == str] #Get final set of introgressions of that strain
  tmpGenes <- tmpMat$Genes
  f <- f[f$SAPAorthID %in% tmpGenes,] #Only retain those introgressions
  no <- length(unique(f$SAPAorthID)) #dim(f)[1]
  nb <- length(unique(tmpMat$Introg_name))
  strainV <- c(strainV,str)
  noIntrosV <- c(noIntrosV,no)
  noBlocksV <- c(noBlocksV,nb)
  f_filtrado <- f[f$Ratio > 0.25, ]
  f_unico <- f_filtrado[!duplicated(f_filtrado$SAPAorthID), ]
  heterozygosityV <- c(heterozygosityV,length(f_unico$Ratio)/no)
}


data <- data.frame(ID = strainV, noBlocks = noBlocksV,noIntros = noIntrosV,heterozygosity = heterozygosityV)
data$meanGperB <- data$noIntros/data$noBlocks
data[is.na(data$meanGperB),]$meanGperB <- 0

smpData <- fread(samples)
data <- merge(data,smpData)
setDT(data)

data[is.na(data$heterozygosity)]$heterozygosity <- 0
data$heterozygosity2 <- data$heterozygosity
data[data$noIntros < 2,]$heterozygosity2 <- 0
data[PhyloGroup_SACE469 !="MMMG" & PhyloGroup_SACE469!="Tamaulipas" & PhyloGroup_SACE469!="Tequila" & PhyloGroup_SACE469!= "SAM2" & PhyloGroup_SACE469!="FG" & PhyloGroup_SACE469!="Alpechin"]$PhyloGroup_SACE469 <- "Others"
wb3 <- c("BR021c1","BR010c1","BR011c1")
data[data$ID %in% wb3,]$PhyloGroup_SACE469 <- "WB3"
data$PhyloGroup_SACE469 <- factor(data$PhyloGroup_SACE469, levels <- c("Others","Alpechin","FG","SAM2","WB3","Tequila","MMMG","Tamaulipas"))

#setNames(data$noIntros,data$ID)

p <- ggplot(data,aes(noIntros,PhyloGroup_SACE469)) +
  geom_violin() + theme_classic() +
  geom_jitter(height=0.2,aes(color=heterozygosity2),size=1.8,alpha=0.7,shape=16) +
  scale_color_gradient(low= "yellow",high= "red3")

ggsave(out,plot=p,width = 9, height = 6)

fwrite(file = hetOut,data[,c("ID","noBlocks","meanGperB","noIntros","heterozygosity","heterozygosity2")])
