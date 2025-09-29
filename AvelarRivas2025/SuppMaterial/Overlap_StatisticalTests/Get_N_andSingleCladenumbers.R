# From Analisis_filtered_final_introMatrix.R

genes_en_Total <- unique(filtered_Intro_df$Gene) 
length(genes_en_Total)

subset_dfma <- filtered_Intro_df[filtered_Intro_df$Phylogenetic_Clade %in% c("Mexican_Agave_1"), ]
subset_dfalp <- filtered_Intro_df[filtered_Intro_df$Phylogenetic_Clade %in% c("Alpechin"), ]
subset_df2 <- filtered_Intro_df[filtered_Intro_df$Phylogenetic_Clade %in% c("Mexican_Agave_2", "WB3", "Tequila_Distillery","French Guiana","SAM2"), ]

# Paso 2: Obtener todos los genes que aparecen en esos clados
genes_en_ma <- unique(subset_dfma$Gene)
genes_en_alp <- unique(subset_dfalp$Gene)
genes_en_NT_noMA1 <- unique(subset_df2$Gene)

# Paso 5: Eliminar esos genes — quedarnos solo con los que están exclusivamente en los clados deseados
genes_ma_alp=intersect(genes_en_ma, genes_en_alp)
genes_exclusivos <- setdiff(genes_ma_alp, genes_en_NT_noMA1)
unique(genes_exclusivos)

# Paso 6: Filtrar el dataframe original para quedarte solo con esas introgresiones
filtered_genes_df <- filtered_Intro_df[filtered_Intro_df$Gene %in% genes_exclusivos & 
                                         filtered_Intro_df$Phylogenetic_Clade %in% c("Mexican_Agave_1", "Alpechin"), ]

# Resultado
dim(filtered_genes_df)
##### 
unsubset <- filtered_Intro_df[filtered_Intro_df$Phylogenetic_Clade %in% c( "Alpechin"), ]
unique(unsubset$Gene)


clados_interes <- c("Mexican_Agave_1", "Mexican_Agave_2", 
                    "Tequila_Distillery", "WB3", 
                    "French Guiana", "SAM2", "Alpechin")

# Filtrar solo esos clados
filtered <- filtered_Intro_df[filtered_Intro_df$Phylogenetic_Clade %in% clados_interes, ]

# Contar genes únicos por clado
library(dplyr)
conteo_genes_por_clado <- filtered %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(NumGenes = n_distinct(Gene)) %>%
  arrange(match(Phylogenetic_Clade, clados_interes))

print(conteo_genes_por_clado)
