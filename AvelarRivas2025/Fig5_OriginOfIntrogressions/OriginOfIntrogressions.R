# From Origin_Introgressions_fromVCF.R
# Toma los "^paired_distances_sapa_PLINK_YPS138_" que vienen de plinkDistMat2.sge
# También incluye los datos de las filogenias

#### MAKE origins_extended_Blocks ######

files3 <- list.files(pattern = "^paired_distances_sapa_PLINK_YPS138_")

# este busca las distancias maximas, es mejor para el IBS
con=0
for (arch in files3) { #
  print(con)
  con=con+1
  if (file.exists(arch) && file.info(arch)$size > 0) {
    a <- read.table(arch, sep = ",")
  } else {
    message("El archivo no existe o está vacío: ", arch)
    next
  }
  blockname=substr(arch,29,61)#36-68
  print(blockname)
  blocksize=Bloques[Bloques$Name==blockname,"SizeNT"]
  #a=read.table("paired_distances_sapa_YPS138_08G00270.csv", sep=",")
  #match( lassaces, unique(c(a$V1, a$V2))  )
  #preestassaces=lassaces[match( lassaces, unique(c(a$V1, a$V2))  )]
  #estassaces = preestassaces[!is.na(preestassaces)]
  
  Columnas=unique(metadata$Genetic_Group_2024)
  filas=lassaces
  this_introgression <- data.frame(matrix(ncol = length(Columnas), nrow = length(filas)))
  colnames(this_introgression) <- Columnas
  rownames(this_introgression) <- filas
  this_introgression_dist <- data.frame(matrix(ncol = length(Columnas), nrow = length(filas)))
  colnames(this_introgression_dist) <- Columnas
  rownames(this_introgression_dist) <- filas
  
  for (i in lassaces) { # lassaces # estassaces
   #i="JS150c1"
    suba = a[a$V1 == i | a$V2 == i, ]
    if (dim(suba)[1] > 0) {
      print(i)
      filtered_data <- suba %>%
        filter((V1 %in% metadata$ID | V2 %in% metadata$ID) & !is.na(V3))
      dim(filtered_data)
      
      # CAMBIO: usamos max en lugar de min
      subamax = filtered_data[filtered_data$V3 == max(filtered_data$V3), ]
      distMax = max(filtered_data$V3)
      
      unicos = unique(c(subamax$V1, subamax$V2))
      losbuenos = match(unicos, metadata$ID)
      metadata$Genetic_Group_2024[losbuenos]
      
      xxx = table(metadata$Genetic_Group_2024[losbuenos])
      xxxdf = as.data.frame.table(xxx)
     
      for (j in 1:dim(xxxdf)[1]) {
        this_introgression[i, as.character(xxxdf$Var1[j])] = as.numeric(as.character(xxxdf$Freq[j]))
        this_introgression_dist[i, as.character(xxxdf$Var1[j])] = distMax / blocksize
      } 
    } 
  }print(xxxdf)
  tamano = nchar(arch)
  write.csv(this_introgression, paste("wblocks/origins_extended_Blocks", substr(arch, 36, nchar(arch)), sep = ""))
  
}



### CARGAR BLOQUES ####
Bloques=read.table("All_blocks_.csv", sep=",",header=TRUE)
Bloques=read.table("C:/Users/HP/Desktop/FromIvan/All_blocks_250609.csv", sep=",",header=TRUE)
Bloques2=read.table("C:/Users/HP/Desktop/FromIvan/All_blocks Oct2024.csv", sep=",",header=TRUE)

Bloques_comb <- rbind(Bloques, Bloques2)
Bloques_comb_unique <- unique(Bloques_comb)
Bloques=Bloques_comb_unique

unlist(Bloques_comb$Genes)
all_genes <- sort(unique(unlist(strsplit(Bloques_comb$Genes, ";"))))

###### BLOQUES aquí va agarrando los origins_extended y genera: Strict_origins NonStrict_origins AllNonStrict_origins ######

files4 <- list.files("wblocks", pattern = "^origins_extended_")
GenesBlocks=substr(all_genes,8,15)
GenesBlocks[1:7]=all_genes[1:7]
filas=lassaces
for (hide in 1) { # generar los dfs donde se van a ir guardando cosas
  Strict_origins_Bl <- data.frame(matrix(ncol = length(GenesBlocks), nrow = length(filas)))
  colnames(Strict_origins_Bl) <- GenesBlocks
  rownames(Strict_origins_Bl) <- filas
  NonStrict_origins_Bl <- data.frame(matrix(ncol = length(GenesBlocks), nrow = length(filas)))
  colnames(NonStrict_origins_Bl) <- GenesBlocks
  rownames(NonStrict_origins_Bl) <- filas
  AllNonStrict_origins_Bl <- data.frame(matrix(ncol = length(GenesBlocks), nrow = length(filas)))
  colnames(AllNonStrict_origins_Bl) <- GenesBlocks # substr(files2,18,25)
  rownames(AllNonStrict_origins_Bl) <- filas
  #arch= "origins_extended_06G00080.csv" #JS933c1  substr(arch, 36, nchar(arch)-4 )
}

con=0
#for(arch in files4){ # recorre todos los archivos
  con=con+1
  print(arch)
  b=read.table(paste("wblocks/",arch,sep=""), sep=",",header=TRUE)
  #b=read.table(arch, sep=",",header=TRUE)
  verdaderos=rowSums(b[,2:17], na.rm=TRUE)>0
  dim(b[verdaderos,])
  subb=b[verdaderos,]
  print(rowSums(subb[,2:17], na.rm = TRUE))
  #gen1=substr(arch, 25,32)
  #gen2=substr(arch, 43,nchar(arch)-4 )
  #gen1=substr(arch, 31,38)
  #gen2=substr(arch, 49,nchar(arch) )
  gen1=substr(arch, 24,31)
  gen2=substr(arch, 42,nchar(arch)-4 )
  if ( nchar(gen2)==0 ){
    gen2=gen1
  }
  subb$SpB_All <- rowSums(subb[, c("SpB", "SpB_1", "YMX_Oax")], na.rm = TRUE)
  subb$SpB_NE <- rowSums(subb[, "YMX_NE", drop = FALSE], na.rm = TRUE) # subb[, "YMX_NE", drop = FALSE]
  subb$SpA_All <- rowSums(subb[, c("SpA", "YMX_SpA")], na.rm = TRUE)
  subb$SpD_All <- rowSums(subb[, c("SpD_1", "SpD")], na.rm = TRUE)
  subb$SpB_Mx <- rowSums(subb[, "YMX_tank", drop=FALSE], na.rm = TRUE)
  subb$SpB_Dgo <- rowSums(subb[, "YMX_Dgo", drop=FALSE], na.rm = TRUE)
  df <- subb[, !colnames(subb) %in% c("SpB", "SpB_1", "YMX_Oax","SpA", "YMX_SpA","SpD_1", "SpD",
                                      "YMX_NE", "YMX_tank", "YMX_Dgo")]
  #if sum(which(colnames(Strict_origins_Bl)==gen1))>0 &  {  }
  rango = which(colnames(Strict_origins_Bl)==gen1) : which(colnames(Strict_origins_Bl)==gen2)
  
  #print(df)
  df[df == 0] <- NA
  if (dim(df)[1]>0) {
    colus=dim(df)[2]
    for ( i in 1:dim(df)[1] ){
      #i=3
      print(i)
      if ( sum(!is.na(df[i,2:colus])) ==1) {
        #nombres=colnames(df)[which(max(df[i,2:colus], na.rm=TRUE)==df[i,])]
        Strict_origins_Bl[df$X[i],rango]=paste( colnames(df)[which(max(df[i,2:colus], na.rm=TRUE)==df[i,])], collapse = " ")
        NonStrict_origins_Bl[df$X[i],rango]=paste( colnames(df)[which(max(df[i,2:colus], na.rm=TRUE)==df[i,])], collapse = " ")
        AllNonStrict_origins_Bl[df$X[i],rango]=paste( colnames(df)[which(max(df[i,2:colus], na.rm=TRUE)==df[i,])], collapse = " ")
      } else {
        NonStrict_origins_Bl[df$X[i],rango]=paste( colnames(df)[which(max(df[i,2:colus], na.rm=TRUE)==df[i,])], collapse = " ")
        AllNonStrict_origins_Bl[df$X[i],rango]=paste( colnames(df)[ 1+which( df[i,2:colus]>0 ) ], collapse = " ")
      }
    }
  }
}

write.csv(Strict_origins_Bl, "Strict_origins_339x1315_Bl_250613.csv" ) #Strict_origins_Blocks
write.csv(AllNonStrict_origins_Bl, "AllNonStrict_origins_339x1315_Bl_250613.csv" )
write.csv(NonStrict_origins_Bl, "NonStrict_origins_339x1315_Bl_250613.csv" )



##### unirlo con merged_genes_expanded (Output del análisis de las filogeinas) #####
library(tidyr)
library(dplyr)

Strict_origins_Bl=read.csv("Strict_origins_339x1315_Bl_250610.csv")
rownames(Strict_origins_Bl)=Strict_origins_Bl$X
Strict_origins_Bl <- Strict_origins_Bl %>% select(-X)
colnames(Strict_origins_Bl) <- sub("^X", "", colnames(Strict_origins_Bl))

Strict_origins_Bl_long <- Strict_origins_Bl %>%
  tibble::rownames_to_column(var = "Strain") %>%  # Convertir nombres de fila en columna
  pivot_longer(-Strain, names_to = "Gene", values_to = "ClosestCladeSimilarity") %>%
  filter(!is.na(ClosestCladeSimilarity))  # Eliminar NA

dim(Strict_origins_Bl_long)
dim(merged_genes_expanded) # este viene de AnalyzeIntrogressionTrees.R

merged_with_similarity <- merged_genes_expanded %>%
  left_join(Strict_origins_Bl_long, by = c("Strain", "Gene"))
View(merged_with_similarity)

filtered_Intro_df <- merged_with_similarity %>%
  filter(TrueIntrogressions == TRUE,
         !is.na(TrueOrigin),
         TrueOrigin != "Sp_und_outgroup",
         TrueOrigin != "WINE-JS445c1SAPA",)
#write.csv(filtered_Intro_df, "StrictOrigins_mergedGenesExpanded_CompleteTable.csv")

table(filtered_Intro_df$Strain)

strain_clade <- merged_genes_expanded %>%
  select(Strain, Phylogenetic_Clade) %>%
  distinct()

# Generar la tabla con los conteos
tabla_por_strain <- filtered_Intro_df %>%
  group_by(Strain, ClosestCladeSimilarity) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = ClosestCladeSimilarity, values_from = n, values_fill = 0) %>%
  left_join(strain_clade, by = "Strain") %>%
  relocate(Phylogenetic_Clade, .after = Strain)

View(tabla_por_strain)

tabla_larga <- tabla_por_strain %>%
  pivot_longer(
    cols = -c(Strain, Phylogenetic_Clade),
    names_to = "ClosestCladeSimilarity",
    values_to = "Count"
  )

stats_df <- tabla_larga %>%
  group_by(Phylogenetic_Clade, ClosestCladeSimilarity) %>%
  summarise(
    Average = mean(Count, na.rm = TRUE),
    StdDev = sd(Count, na.rm = TRUE),
    .groups = "drop"
  )

facet_order <- c("Mexican_Agave_1", "Mexican_Agave_2", "Tequila_Distillery",  "WB3","SAM2","French Guiana", "Alpechin")
stats_df$Phylogenetic_Clade <- factor(stats_df$Phylogenetic_Clade, levels = facet_order)
ordered_origins <- c("SpB_Mx", "SpB_Dgo", "SpB_NE", "SpB_All", "SpB_Brasil",
                     "SpD_All", "SpC", "SpC.", "Hawaii", "SpA_All", "FarEastAsia","NA")

# Convertir la columna ClosestCladeSimilarity en factor con orden
stats_df$ClosestCladeSimilarity <- factor(stats_df$ClosestCladeSimilarity, levels = ordered_origins)

# PLOT de NUMERO DE ORIGENES POR CLADO
stats_df_filtered <- stats_df %>%
  filter(Phylogenetic_Clade != "Alpechin" & !is.na(Phylogenetic_Clade))

pdf("Origin_perClade_free_y_scale3.pdf",5,5)
ggplot(stats_df_filtered , aes(x = ClosestCladeSimilarity, y = Average)) +
  geom_bar(stat = "identity", fill = "grey50", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = Average - StdDev, ymax = Average + StdDev),
                position = position_dodge(width = 0.7), width = 0.2) +
  facet_wrap(~ Phylogenetic_Clade, ncol = 1, scales = "free_y") +
  labs(x = "Origin", y = "Introgressed genes per strain ± SD") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


##### For the final revision  ########
library(dplyr)
library(ggplot2)

facet_order_Names <- c("MA1", "MA2", "TeqD", "SAM2", "FG", "WB3")
names(facet_order_Names) <- c("Mexican_Agave_1",
                              "Mexican_Agave_2",
                              "Tequila_Distillery",
                              "SAM2",
                              "French Guiana",
                              "WB3")

ordered_origins       <- c("SpB_Mx", "SpB_Dgo", "SpB_NE", "SpB_All", "SpB_Brasil")
ordered_origins_Names <- c("SpB_MxAg", "SpB_Mx2", "SpB_Mx1", "SpB", "SpB_Bra")
names(ordered_origins_Names) <- ordered_origins


# Datos crudos filtrados (sin NA y solo esos 5 origins)
tabla_larga_filtrada <- tabla_larga %>%
  filter(
    Phylogenetic_Clade %in% names(facet_order_Names),  # <-- solo los clados que quieres
    !is.na(Phylogenetic_Clade),
    !is.na(ClosestCladeSimilarity),
    ClosestCladeSimilarity %in% ordered_origins
  ) %>%
  mutate(
    Phylogenetic_Clade     = factor(Phylogenetic_Clade, levels = facet_order),
    ClosestCladeSimilarity = factor(ClosestCladeSimilarity,
                                    levels = rev(ordered_origins))
  ) %>%
  droplevels()

# Resumen (media, SD y N)
stats_df <- tabla_larga_filtrada %>%
  group_by(ClosestCladeSimilarity, Phylogenetic_Clade) %>%
  summarise(
    Average = mean(Count, na.rm = TRUE),
    StdDev  = sd(Count, na.rm = TRUE),
    N       = sum(!is.na(Count)),
    .groups = "drop"
  ) %>%
  droplevels()

# n por clado para los títulos de los facets
n_by_clade <- tabla_larga_filtrada %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(N_clade = n_distinct(Strain), .groups = "drop")

# nombres cortos para los facets
facet_order_Names <- c("MA1", "MA2", "TeqD", "SAM2", "FG", "WB3")
names(facet_order_Names) <- c("Mexican_Agave_1", "Mexican_Agave_2", "Tequila_Distillery", "SAM2", "French Guiana", "WB3")

facet_labels <- setNames(
  paste0(facet_order_Names[n_by_clade$Phylogenetic_Clade], " (n=", n_by_clade$N_clade, ")"),
  n_by_clade$Phylogenetic_Clade
)

ordered_origins_Names <- c("SpB_MxAg", "SpB_Mx2", "SpB_Mx1", "SpB", "SpB_Bra")
names(ordered_origins_Names) <- ordered_origins

max_xlim <- 280

pdf("Fig5a_BarsAndPlots_b.pdf", 88/25.4, 80/25.4 )
ggplot() +
  geom_col(data = stats_df, aes(x = Average, y = ClosestCladeSimilarity),
           fill  = "grey50", width = 0.7) +
  geom_errorbar(data = stats_df,
                aes(xmin = Average - StdDev, xmax = Average + StdDev, y = ClosestCladeSimilarity),
                width = 0.2) +
  geom_jitter(data = tabla_larga_filtrada,
              aes(x = Count, y = ClosestCladeSimilarity),
              height = 0.15, alpha  = 0.2, size   = 1) +
  facet_wrap(~ Phylogenetic_Clade,
             ncol = 7, scales = "fixed", labeller = as_labeller(facet_labels)) +
  scale_y_discrete(labels = ordered_origins_Names) +
  coord_cartesian(xlim = c(0, max_xlim)) +
  labs(x = "Introgressed genes per strain (mean ± SD)",
       y = "Origin") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_text(size = 8))
dev.off()

