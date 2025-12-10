################################################################################
##############                  LIBRARIES NEEDED                  ##############
################################################################################
library(data.table)
library(ggh4x)
library(ggvenn)
################################################################################
##############                       INPUTS                       ##############
################################################################################
peterMatrixFile <- "./IntrogressedGenesProcessed.tsv")
ygmxMatrixFile <- "../../figure4/fig4A/pipeline/trueIntrogressions/data/ascTrees/updatedMatrix020725_final.csv"
telliniMatrixFile <- "./TelliniIntrogressedGenes.csv"
equivalenceTelliniFile <- "./equivalenceTellini.tsv"
equivalencePeterFile <- "./equivalencePeter.tsv"
equivalenceYGMXFile <- "./equivalenceYGMX.tsv"
YGMXGenesDictFile <- "./GenesEquivalenceYPS138_v2.txt"
runIdsFile <- "./idNamesRelation.csv"
commonstrainsFile <- "./strainsInCommon.txt"
################################################################################
##############                      OUTPUTS                       ##############
################################################################################
outPng="./compareTelliniPeter.pdf"
################################################################################
##############                     MAIN CODE                      ##############
################################################################################
#Load  file with introgressions
##Provided by Mateo from Josepsh lab
peterMatrix <- fread(peterMatrixFile)
##Obtaine for figure 3C
ygmxMatrix <- fread(ygmxMatrixFile)
#Change column name
#colnames(ygmxMatrix)[2] <- "Phylogenetic_Clade"
#From Abraham, internal control
runIds <- read.csv(runIdsFile)
#Change from Ids to Sample Names
ygmxMatrix <- merge(runIds[,c(1,4)],ygmxMatrix)
ygmxMatrix$ID <- NULL
#Erase incorrectly asigned strains
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "AID",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "AIE",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "AHM",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "ABM",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "AQA",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "AQD",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "AVB",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "ATQ",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "BRK",]
ygmxMatrix <- ygmxMatrix[!ygmxMatrix$Sample_Name == "BFN",]
##Supp. table 6 from tellini 2024
telliniMatrix <- fread(telliniMatrixFile)
#Rename some columns that arise problems
colnames(telliniMatrix)[colnames(telliniMatrix) == 'X11_3002'] = '11_3002'
colnames(telliniMatrix)[colnames(telliniMatrix) == 'X5_3004'] = '5_3004'
colnames(telliniMatrix)[colnames(telliniMatrix) == 'X59A'] = '59A'
colnames(telliniMatrix)[colnames(telliniMatrix) == 'X6320_A7'] = '6320_A7'
colnames(telliniMatrix)[colnames(telliniMatrix) == 'X2_1392'] = '2_1392'
colnames(telliniMatrix)[colnames(telliniMatrix) == 'X229'] = '229'
#Load files with name equivalences
##From supplementary table 2 from Tellini et al 2024
equivalenceTellini <- fread(equivalenceTelliniFile,head=FALSE,col.names=c("Universal","Tellini"))
##From supplementary table 1 from peter et al 2018
equivalencePeter <- fread(equivalencePeterFile,head=FALSE,col.names=c("Universal","Peter"))
##From supplementary table 2 from avelar-rivas et al 2025
equivalenceYGMX <- fread(equivalenceYGMXFile,head=FALSE,col.names=c("Universal","Ygmx"))

ygmxDict <- fread(YGMXGenesDictFile,col.names=c("ygmx","Peter"))
##Change "-" in names of strains for "_". This causes trouble
equivalenceTellini$Universal <- gsub("-","_",equivalenceTellini$Universal)
equivalencePeter$Universal <- gsub("-","_",equivalencePeter$Universal)
equivalenceYGMX$Universal <- gsub("-","_",equivalenceYGMX$Universal)
#From the ygmx genes dictionary we will retain only those with a non dubious sgd ortholog.
ygmxDictFiltered <- ygmxDict[Peter != "NA"]
ygmxDictFiltered <-ygmxDictFiltered[!grepl("/", ygmxDictFiltered$Peter),]
#From the peter genes we will retain only those that have an annotation.
peterFiltered <- peterMatrix[`Ortholog in SGD_2010` != ""][,-c(1)]
#Delete row with NA's, should talk with mateo dechiara
peterFiltered <- peterFiltered[complete.cases(peterFiltered),]
#From the peter genes we will retain only those with a nondubious sgd ortholog.
peterFiltered <- peterFiltered[!grepl(",", peterFiltered$`Ortholog in SGD_2010`)]
#ygmx genes matrix is transform so instead of origins we have presence and absence.
#Presence is encoded as integer for uniformity with the other data.
TygmxMatrix <- as.data.frame(t(ygmxMatrix[,2:ncol(ygmxMatrix)]),row.names=FALSE)
#TygmxMatrix[!is.na(TygmxMatrix)] <- 1
#TygmxMatrix[is.na(TygmxMatrix)] <- 0
#TygmxMatrix[] <- lapply(TygmxMatrix, as.integer)
TygmxMatrix$`ygmxTmp` <- colnames(ygmxMatrix)[-c(1)]
#TygmxMatrix[!grepl("YPS", TygmxMatrix$`ygmxTmp`),]$`ygmxTmp` <- paste0("YPS_",TygmxMatrix[!grepl("YPS", TygmxMatrix$`ygmxTmp`),]$`ygmxTmp`)
namesListed <- strsplit(TygmxMatrix$ygmxTmp,"_")
TygmxMatrix$`ygmx` <- unlist(namesListed)[2*(1:length(namesListed))]
TygmxMatrix$`ygmxTmp` <- NULL
colnames(TygmxMatrix) <- c(ygmxMatrix$Sample_Name,"ygmx")
#We filter our data, retaining only genes with non dubious annotations.
ygmxFiltered <- merge(TygmxMatrix,ygmxDictFiltered)[-c(1)]
#The genes from tellini are transformed in a way such that only genes with 60% or more
#are considered present
telliniCutoff <- telliniMatrix[,-1]
telliniCutoff[telliniCutoff < 60] <- 0
telliniCutoff[telliniCutoff >= 60] <- 1
telliniCutoff$Peter <- telliniMatrix[,1]
#Data frames with equivaleces of names are constructed in a way such that all data can
#be compared.
namesTellini <- data.frame(colnames(telliniCutoff[,-c("Peter")]))
colnames(namesTellini) <- "Tellini"
equivNamesTell <- merge(namesTellini,equivalenceTellini)
namesPeter <- data.frame(colnames(peterFiltered[,-c("Ortholog in SGD_2010")]))
colnames(namesPeter) <- "Peter"
equivNamesPeter <- merge(namesPeter,equivalencePeter)
namesYgmx <- data.frame(colnames(ygmxFiltered[,-ncol(ygmxFiltered)]))
colnames(namesYgmx) <- "Ygmx"
equivNamesYgmx <- merge(namesYgmx,equivalenceYGMX)

strainsAll <- read.csv(commonstrainsFile,header = FALSE)$V1

#Reorder columns alphabetically
telliniColnames <- colnames(telliniCutoff[,-c("Peter")])
telliniCutoff <- telliniCutoff[,c(telliniColnames[order(telliniColnames)],"Peter"),with = FALSE]
peterColnames <- colnames(peterFiltered[,-c("Ortholog in SGD_2010")])
peterFiltered <- peterFiltered[,c("Ortholog in SGD_2010",peterColnames[order(peterColnames)]),with = FALSE]
ygmxColnames <- colnames(ygmxFiltered[,-ncol(ygmxFiltered)])
ygmxFiltered <- ygmxFiltered[,c(ygmxColnames[order(ygmxColnames)],"Peter")]
#Change name of columns to make them congruent
orderedTellini <- equivNamesTell$Universal[order(equivNamesTell$Tellini)]
orderedPeter <- equivNamesPeter$Universal[order(equivNamesPeter$Peter)]
orderedYgmx <- equivNamesYgmx$Universal[order(equivNamesYgmx$Ygmx)]
colnames(telliniCutoff) <- c(orderedTellini,"Gene")
colnames(peterFiltered) <- c("Gene",orderedPeter)
colnames(ygmxFiltered) <- c(orderedYgmx,"Gene")
#Delete all duplicates
telliniCutoff <- telliniCutoff[!(duplicated(telliniCutoff$Gene) | duplicated(telliniCutoff$Gene, fromLast = TRUE)),]
peterFiltered <- peterFiltered[!(duplicated(peterFiltered$Gene) | duplicated(peterFiltered$Gene, fromLast = TRUE)),]
ygmxFiltered <- ygmxFiltered[!(duplicated(ygmxFiltered$Gene) | duplicated(ygmxFiltered$Gene, fromLast = TRUE)),]
#Recover only the strains shared by all strains
telliniUseful <- telliniCutoff[,c("Gene",strainsAll),with=FALSE]
peterUseful <- peterFiltered[,c("Gene",strainsAll),with=FALSE]
ygmxUseful <- data.table(ygmxFiltered[,c("Gene",strainsAll)])

viy <- c()
vip <- c()
vit <- c()
genevector <- unique(intersect(intersect(peterUseful$Gene,telliniUseful$Gene),ygmxUseful$Gene))
for(strain in strainsAll){
  for(gene in ygmxUseful$Gene){
    if(ygmxUseful[ygmxUseful$Gene == gene,][[strain]]){
      viy <- c(viy,paste(gene,strain,sep = "@"))
    }
  }
  for(gene in peterUseful$Gene){
    if(peterUseful[peterUseful$Gene == gene,][[strain]]){
      vip <- c(vip,paste(gene,strain,sep = "@"))
    }
  }
  for(gene in telliniUseful$Gene){
    if(telliniUseful[telliniUseful$Gene == gene,][[strain]]){
      vit <- c(vit,paste(gene,strain,sep = "@"))
    }
  }
}
yptI = 0
for(strain in strainsAll){
  for(gene in genevector){
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    if(ygmxUseful[ygmxUseful$Gene == gene,][[strain]]){tmp1 = 1}
    if(peterUseful[peterUseful$Gene == gene,][[strain]]){tmp2 = 1}
    if(telliniUseful[telliniUseful$Gene == gene,][[strain]]){tmp3 = 1}
    if(tmp1 == 1 & tmp2 == 1 & tmp3 == 1){yptI = yptI + 1}
  }
}
vectorIntrogsYgmx <- viy
vectorIntrogsPeter <- vip
vectorIntrogsTellini <- vit
#Make venn diagram
allIntrogressionsVenn <- ggvenn(list(Tellini = vectorIntrogsTellini, Peter = vectorIntrogsPeter, Ygmx = vectorIntrogsYgmx),fill_color = c("red","green","blue"))
#Organize genes present in each strain in lists of vectors
presentTellini <- list()
for(strain in strainsAll){
  tmp <- ifelse(telliniUseful[,..strain] == 1,telliniUseful$Gene,0)
  tmp <- tmp[!tmp==0]
  presentTellini[[strain]][["Genes"]] <- tmp
  presentTellini[[strain]][["Number"]] <- length(tmp)
}
presentPeter <- list()
for(strain in strainsAll){
  tmp <- ifelse(peterUseful[,..strain] == 1,peterUseful$Gene,0)
  tmp <- tmp[!tmp==0]
  presentPeter[[strain]][["Genes"]] <- tmp
  presentPeter[[strain]][["Number"]] <- length(tmp)
}
presentYgmx <- list()
for(strain in strainsAll){
  tmp <- ifelse(ygmxUseful[,..strain] == 1,ygmxUseful$Gene,0)
  tmp <- tmp[!tmp==0]
  presentYgmx[[strain]][["Genes"]] <- tmp
  presentYgmx[[strain]][["Number"]] <- length(tmp)
}

allIntrogressionsVenn <- ggvenn(list(Tellini = vectorIntrogsTellini, Peter = vectorIntrogsPeter, Ygmx = vectorIntrogsYgmx),fill_color = c("red","green","blue"))

ggsave(filename=outPng, plot = allIntrogressionsVenn)
