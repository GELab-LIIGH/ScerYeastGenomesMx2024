################################################################################
##############                  LIBRARIES NEEDED                  ##############
################################################################################
library(data.table)
library(ggvenn)
################################################################################
##############                       INPUTS                       ##############
################################################################################
barbosaIntros <- "./barbosaIntros.csv"
ygmxIntros <- "../../figure4/fig4A/pipeline/trueIntrogressions/data/ascTrees/updatedMatrix100425_v2.csv"
telliniIntros <- "./telliniIntros.csv"
equivalenceTelliniFile <- "./equivalenceTellini.tsv"
equivalenceYGMXFile <- "./equivalenceYGMX.tsv"

YGMXGenesDictFile <- "./GenesEquivalenceYPS138_v2.txt"
runIdsFile <- "./idNamesRelation.csv"
################################################################################
##############                      OUTPUTS                       ##############
################################################################################
outPng <- "./compareTelliniBarbosa.pdf"
################################################################################

################################################################################
##############                     MAIN CODE                      ##############
################################################################################
#From barbosa supplementary table
barbosaMatrix <- read.csv(barbosaIntros)
colnames(barbosaMatrix) <- gsub("\\.","-",colnames(barbosaMatrix))
#From our data
ygmxMatrix <- read.csv(ygmxIntros)
#From tellini supplementary data
telliniMatrix <- read.csv(telliniIntros)
colnames(telliniMatrix) <- gsub("\\.","-",colnames(telliniMatrix))
#
runIds <- read.csv(runIdsFile)
#Merge to match IDs
ygmxMatrix <- merge(runIds[,c(1,4)],ygmxMatrix)
ygmxMatrix <- ygmxMatrix[ygmxMatrix$Sample_Name == "UFMG-CM-Y642" | ygmxMatrix$Sample_Name == "UFMG-CM-Y263" | ygmxMatrix$Sample_Name == "UFMG-CM-Y641",]
ygmxMatrix$ID <- NULL
#ygmx genes matrix is transform so instead of origins we have presence and absence.
#Presence is encoded as integer for uniformity with the other data.
TygmxMatrix <- as.data.frame(t(ygmxMatrix[,2:ncol(ygmxMatrix)]),row.names=FALSE)
TygmxMatrix$`ygmxTmp` <- colnames(ygmxMatrix)[-c(1)]
TygmxMatrix <- TygmxMatrix[grep("YPS138",TygmxMatrix$ygmxTmp),]
namesListed <- strsplit(TygmxMatrix$ygmxTmp,"_")
TygmxMatrix$`ygmx` <- unlist(namesListed)[2*(1:length(namesListed))]
TygmxMatrix$`ygmxTmp` <- NULL
colnames(TygmxMatrix) <- c(ygmxMatrix$Sample_Name,"ygmx")

ygmxDict <- fread(YGMXGenesDictFile,col.names=c("ygmx","Peter"))
#From the ygmx genes dictionary we will retain only those with a non dubious sgd ortholog.
ygmxDictFiltered <- ygmxDict[Peter != "NA"]
ygmxDictFiltered <-ygmxDictFiltered[!grepl("/", ygmxDictFiltered$Peter),]

#We filter our data, retaining only genes with non dubious annotations.
ygmxFiltered <- merge(TygmxMatrix,ygmxDictFiltered)[-c(1)]
#Delete all duplicates
ygmxFiltered <- ygmxFiltered[!(duplicated(ygmxFiltered$Peter) | duplicated(ygmxFiltered$Peter, fromLast = TRUE)),]
vectorBarbosa <- c()
vectorTellini <- c()
vectorYgmx <- c()

for(strain in colnames(barbosaMatrix)){
  tmp <- paste0(barbosaMatrix[,strain],"@",strain)
  vectorBarbosa <- unique(c(vectorBarbosa,tmp))
}

for(strain in colnames(telliniMatrix)){
  tmp <- paste0(telliniMatrix[,strain],"@",strain)
  vectorTellini <- unique(c(vectorTellini,tmp))
}

for(strain in ygmxMatrix$Sample_Name){
  tmp <- c()
  for(gene in ygmxFiltered$Peter){
    if(ygmxFiltered[ygmxFiltered$Peter == gene ,strain] == 1){
      tmp <- c(tmp,paste0(gene,"@",strain))}
  }
  vectorYgmx <- unique(c(vectorYgmx,tmp))
}
allIntrogressionsVenn <- ggvenn(list(Tellini = vectorTellini, Barbosa = vectorBarbosa, Ygmx = vectorYgmx),fill_color = c("red","green","blue"))


ggsave(filename=outPng, plot = allIntrogressionsVenn)
