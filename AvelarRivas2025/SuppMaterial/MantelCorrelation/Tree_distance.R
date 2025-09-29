install.packages("ape")
install.packages("dplyr")
install.packages("tidyr")
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)

tree<-read.tree("../RAxML_bipartitionsBranchLabels.Matrix_SNPs_SACE_from_CONC_gt_onlySNPs_filtered_missing_10_biallelic_2.tree")
tree2<-read.tree("C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/Supp_Material/BuildSuppTable/tree_SACE467_TreeName.nwk")

PatristicDistMatrix<-cophenetic(tree2)
dim(PatristicDistMatrix)

new_Mgd <- read.csv("D:/Dropbox/Posdoc/lavis/int/AllSegments/AllIntrogresions/merged_df_ISandGroups3.csv")

fg=match( new_Mgd[new_Mgd$PhyloGroup_SACE469=="FG","Row_Names"], colnames(PatristicDistMatrix) )
FG=fg[!is.na(fg)]

mmmg=match( new_Mgd[new_Mgd$PhyloGroup_SACE469=="MMMG","Row_Names"], colnames(PatristicDistMatrix) )
MA2=mmmg[!is.na(mmmg)]

tmps=match( new_Mgd[new_Mgd$PhyloGroup_SACE469=="Tamaulipas","Row_Names"], colnames(PatristicDistMatrix) )
MA1=tmps[!is.na(tmps)]


MA1_MA2=PatristicDistMatrix[MA1,MA2]
FG_MA2=PatristicDistMatrix[FG,MA2]
FG_MA1=PatristicDistMatrix[FG,MA1]
MA1_MA1=PatristicDistMatrix[MA1,MA1]
MA2_MA2=PatristicDistMatrix[MA2,MA2]
FG_FG=PatristicDistMatrix[FG,FG]


aMA1_MA2 <- data.frame(value = as.vector(as.matrix(MA1_MA2)), dataframe = "MA1_MA2")
bFG_MA2 <- data.frame(value = as.vector(as.matrix(FG_MA2)), dataframe = "FG_MA2")
cFG_MA1 <- data.frame(value = as.vector(as.matrix(FG_MA1)), dataframe = "FG_MA1")

aMA1_MA1 <- data.frame(value = as.vector(as.matrix(MA1_MA1)), dataframe = "MA1_MA1")
bMA2_MA2 <- data.frame(value = as.vector(as.matrix(MA2_MA2)), dataframe = "MA2_MA2")
cFG_FG <- data.frame(value = as.vector(as.matrix(FG_FG)), dataframe = "FG_FG")




combined_df <- bind_rows(aMA1_MA2, bFG_MA2, cFG_MA1, aMA1_MA1, bMA2_MA2, cFG_FG)


ggplot(combined_df, aes(x = dataframe, y = value)) +
  geom_violin() +
  theme_minimal() +
  labs(x = "DataFrame",
       y = "Valor")
###  CALCULAR MATRIZ DE DISTANCIAS  ####
install.packages("geosphere")
library(geosphere)
suptable=read.csv("C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/Supp_Material/Supplementary_Table1.csv")
coord1<-c(suptable$Longitude[1],suptable$Latitude[1])
coord2<-c(suptable$Longitude[10],suptable$Latitude[10])


submat=suptable[(!is.na(suptable$Latitude)),c("Sample_Name","Latitude","Longitude")]
coordenadas<-as.matrix(submat[, c("Longitude", "Latitude")])
distancia <- distHaversine(coordenadas)

submat=suptable[(!is.na(suptable$Latitude)),c("Sample_Name","Latitude","Longitude")]
coordenadas<-as.matrix(submat[, c("Longitude", "Latitude")])
matriz_distancias <- distm(coordenadas, fun = distHaversine)
rownames(matriz_distancias) <- submat$Sample_Name
colnames(matriz_distancias) <- submat$Sample_Name



###   #####
install.packages("vegan")
library(vegan)
sum(!is.na(match(rownames(PatristicDistMatrix), submat$Sample_Name)))

a=suptable$PhylogeneticGroup.SuppFig1.=="Mexican_Agave_1"
b=suptable$PhylogeneticGroup.SuppFig1.=="Mexican_Agave_2"
interesantes=suptable[a|b,"Sample_Name"]
  
common_rownames <- intersect(rownames(matriz_distancias), rownames(PatristicDistMatrix))
comm = intersect(common_rownames, interesantes )
matriz_distancias_common <- matriz_distancias[comm, comm]
PatristicDistMatrix_common <- PatristicDistMatrix[comm, comm]
mt1=mantel.test(matriz_distancias_common, PatristicDistMatrix_common, by = "pairwise", permutations = 9999)
mt1
mantel_pearson_MA1MA2 <- mantel(matriz_distancias_common, PatristicDistMatrix_common, method = "pearson", permutations = 9999)
mantel_spearman_MA1MA2 <- mantel(matriz_distancias_common, PatristicDistMatrix_common, method = "spearman", permutations = 9999)
plot(matriz_distancias_common, PatristicDistMatrix_common)


df <- data.frame(
  matriz_distancias = as.vector(matriz_distancias_common),
  patristic_distancias = as.vector(PatristicDistMatrix_common)
)
pearson_stat <- mantel_pearson_MA1MA2$statistic
pearson_pval <- mantel_pearson_MA1MA2$signif
spearman_stat <- mantel_spearman_MA1MA2$statistic
spearman_pval <- mantel_spearman_MA1MA2$signif

# Crea el plot con ggplot2
p <- ggplot(df, aes(x = matriz_distancias, y = patristic_distancias)) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Agrega la línea de regresión lineal
  labs(
    x = "Matriz de Distancias Geográficas",
    y = "Distancia Patristica",
    title = "Mantel Test Results"
  ) +
  theme_minimal() +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Pearson Mantel: r = ", round(pearson_stat, 3), ", p = ", format(pearson_pval, scientific = TRUE, digits = 2), 
                   "\nSpearman Mantel: r = ", round(spearman_stat, 3), ", p = ", format(spearman_pval, scientific = TRUE, digits = 2)),
    hjust = 1.1, vjust = 1.1, size = 5, color = "black"
  )


#####
print("MA1")
mantel_pearson_MA1$signif
mantel_pearson_MA1$statistic
mantel_spearman_MA1$signif
mantel_spearman_MA1$statistic
print("MA2")
mantel_pearson_MA2$signif
mantel_pearson_MA2$statistic
mantel_spearman_MA2$signif
mantel_spearman_MA2$statistic
print("MA1MA2")
mantel_pearson_MA1MA2$signif
mantel_pearson_MA1MA2$statistic
mantel_spearman_MA1MA2$signif
mantel_spearman_MA1MA2$statistic
