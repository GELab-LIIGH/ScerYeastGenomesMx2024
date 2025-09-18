options(stringsAsFactors = F)
library(reshape2)
################################################################################
###########                            INS                           ###########
ori <- "/Users/invitado/Documents/Estancia5/"
directory <- "/mnt/Timina/lmorales/Public/reducedSpaTrees/pipelineBra2/data/ascTrees/results/"
compareDirectory <- "/mnt/Timina/lmorales/isedeno/IntrogressionComparison/"
dir <- paste0(ori,directory)
compareDir <- paste0(ori,compareDirectory)
annotFile <- paste0(ori,"/mnt/Timina/lmorales/isedeno/articleFigures/final/figure3/fig3A/pipeline/data/SAPA_YPS138_v1.txt")
blockListFile <- paste0(ori,"/mnt/Timina/lmorales/Public/reducedSpaTrees/pipelineBra2/data/introgression_coordinates.txt")
################################################################################
################################################################################
###########                           OUTS                           ###########
outFile <- paste0(dir,"updatedMatrix020725_final.csv")
IntrosFile <- paste0(compareDir,"AllIntros020725_final.csv")
################################################################################

################################################################################
###########            STEP 1: Obtain genes in each block            ###########
annot <- read.csv(annotFile,col.names = c("Chr", "Start", "End", "Feat", "Name", "Note"), sep = "\t")
blockList <- read.csv(blockListFile,sep="\t")

genesPerBlock <- list()
for(row in 1:nrow(blockList)){
  block <- blockList[row,]
  blockName <- block$ID
  chr <- block$Chromosome
  s <- block$Start
  e <- block$End
  genesPerBlock[[blockName]] <- annot[annot$Chr == chr & annot$Start >= s & annot$End <= e,]$Name
}
################################################################################

setwd(dir)

filesList <- list.files(path = ".", pattern = "_results_250409.csv")
length(filesList)
strain_l <- c()
intro_l <- c()
value_l <- c()
empty <- c()
nottrue <- 0
out <- 0
for(file in filesList){
  df <- read.csv(file)
  bf <- dim(df)[1]
  df <- df[df$TrueIntrogressions,]
  af <- dim(df)[1]
  tq <- bf - af
  nottrue <- nottrue + tq
  df <- df[df$TrueOrigin != "Sp_und_outgroup",] #For v2, eliminate
  df <- df[df$TrueOrigin != "WINE-JS445c1SAPA",]
  af2 <- dim(df)[1]
  oq <- af - af2
  out <- out + oq
  strain <- strsplit(file,split="_")[[1]][1]
  if(dim(df)[1] > 1){
    for(nrow in 1:nrow(df)){
      row <- df[nrow,]
      location <- row$File
      blockName <- strsplit(strsplit(location,split = "Matrix_")[[1]][2],split= "_filtered")[[1]][1]
      genes <- genesPerBlock[[blockName]]
      l <- length(genes)
      strain_l <- c(strain_l,rep(strain,l))
      value_l <- c(value_l,rep(1,l))
      intro_l <- c(intro_l,genes)
    }
  }
  else{
    empty <- c(empty,strain)
  }
}
strain_l <- c(strain_l,empty)
l0 <- length(empty)
g0 <- intro_l[1]
intro_l <- c(intro_l,rep(g0,l0))
value_l <- c(value_l,rep(0,l0))

df_intros <- data.frame(ID = strain_l, Intro = intro_l, value = value_l)
intro_mat <- dcast(data = df_intros, formula  = ID~Intro, value.var = "value")
intro_mat[is.na(intro_mat)] <- 0
write.csv(intro_mat,file = outFile,row.names=FALSE)
genes <- colnames(intro_mat[,-c(1)])
genes <- sub("^YPS138_", "", genes)
write.table(genes,file = IntrosFile,row.names=FALSE,col.names=FALSE,quote=FALSE)

setwd("/Users/invitado/Documents/Estancia5/mnt/Timina/lmorales/isedeno/IntrogressionComparison")
df_intros <- df_intros[df_intros$value == 1,]
keepp <- c("ERR1309514","ERR1309020","BR021c1","BR010c1","BR011c1","SRR800827","JS440c1")
df_intros[df_intros$ID %in% keepp,]
