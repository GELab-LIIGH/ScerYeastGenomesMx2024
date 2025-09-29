# From Origin_Introgressions_fromVCF.R
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

# #### From Analysis_Filtered_final_introMatrix.R  
setwd("C:/Users/jabra/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/figs_code/figsR/Paired_Distances_250612")
load("250716_mi_sesion.RData")

IntroHet = read.csv("HetTable_2_240903_conAlp_mod.csv")
View(final_combined_table_2)
library(dplyr)

filtered_Intro_df_with_hetinfo <- filtered_Intro_df %>%
  left_join(
    IntroHet %>%
      select(Strain, SAPA_Gene, SAPA_cov, SAPA_covnorm, SACE_Gene,
             SACE_cov, SACE_covnorm, PCT_Identity, Het, HasOrtholog, Ratio1),
    by = c("Strain" = "Strain", "Genes" = "SAPA_Gene")
  )

filtered_Intro_df_with_hetinfo <- filtered_Intro_df_with_hetinfo %>%
  mutate(Ratio1_category = case_when(
    is.na(Ratio1) ~ "0",
    Ratio1 < 0.25 ~ "1",
    Ratio1 >= 0.25 & Ratio1 <= 4 ~ "2",
    Ratio1 > 4 ~ "3"
  ))


# #### From AnalyzeIntrogressionTrees.R
library(dplyr)
library(tidyr)
library(pheatmap)
install.packages("forcats")

orden_cepas_MA1_combined <- c(
  "JS421c1", "JS863c1", "JS462c1", "JS806c1", "JS880c1", "XA121c6", "XA121c18", 
  "YMX507A06", "YMX506F12", "YMX507A04", "YMX507A05", "YMX507A07", 
  "JS273c1", "JS240c1", "YMX005588", "YMX005566", "YMX005587", "XA126c7", 
  "YMX005565", "XA121c20"
)

# 20 cepas aleatorias de French Guiana
french_guiana_strains <- filtered_Intro_df_with_hetinfo %>%
  filter(Phylogenetic_Clade == "French Guiana") %>%
  distinct(Strain) %>%
  slice_sample(n = 20) %>%
  pull(Strain)

# 20 cepas aleatorias de Mexican_Agave_2
mexican_agave2_strains <- filtered_Intro_df_with_hetinfo %>%
  filter(Phylogenetic_Clade == "Mexican_Agave_2") %>%
  distinct(Strain) %>%
  slice_sample(n = 20) %>%
  pull(Strain)

# Todas las cepas de Tequila_Distillery, SAM2 y WB3
tequila_strains <- filtered_Intro_df_with_hetinfo %>%
  filter(Phylogenetic_Clade == "Tequila_Distillery") %>%
  distinct(Strain) %>%
  pull(Strain)

sam2_strains <- filtered_Intro_df_with_hetinfo %>%
  filter(Phylogenetic_Clade == "SAM2") %>%
  distinct(Strain) %>%
  pull(Strain)

wb3_strains <- filtered_Intro_df_with_hetinfo %>%
  filter(Phylogenetic_Clade == "WB3") %>%
  distinct(Strain) %>%
  pull(Strain)

# Orden total de cepas
all_strains_order <- c(
  orden_cepas_MA1_combined,
  mexican_agave2_strains,
  tequila_strains,
  sam2_strains,
  wb3_strains,
  french_guiana_strains#,
  #Alp_strains
)

# ----------------------
# Preparar orden de genes por Chr y Start
# ----------------------

# Extraer orden de genes
gene_order_info <- filtered_Intro_df_with_hetinfo %>%
  filter(!is.na(Genes), !is.na(Chr), !is.na(Start)) %>%
  distinct(Genes, Chr, Start) %>%
  mutate(
    Chr_num = as.numeric(gsub("SAPA_YPS138_v1_chr_", "", Chr))
  ) %>%
  arrange(Chr_num, Start)

# Lista ordenada de genes
ordered_genes <- gene_order_info$Genes

# Crear anotación para cromosoma
annotation_col <- gene_order_info %>%
  select(Genes, Chr_num) %>%
  distinct() %>%
  column_to_rownames("Genes")

# Asegurar factor con niveles del 1 al 16
annotation_col$Chr_num <- factor(annotation_col$Chr_num, levels = 1:16)

# ----------------------
# Crear matriz de Ratio1
# ----------------------
library(tibble)
gene_matrix_all <- filtered_Intro_df_with_hetinfo %>%
  filter(Strain %in% all_strains_order) %>%
  filter(!is.na(Genes), !is.na(Ratio1)) %>%
  group_by(Strain, Genes) %>%
  summarise(Ratio1 = max(Ratio1), .groups = "drop") %>%
  pivot_wider(
    names_from = Genes,
    values_from = Ratio1,
    values_fill = NA
  ) %>%
  slice(match(all_strains_order, Strain)) %>%
  column_to_rownames("Strain")

# Asegurar que las columnas estén en el orden correcto
gene_matrix_all <- gene_matrix_all[, intersect(ordered_genes, colnames(gene_matrix_all))]

# Reemplazar NA por -0.1
gene_matrix_all[is.na(gene_matrix_all)] <- -0.1

# ----------------------
# Heatmap
# ----------------------

# Colores y breaks
breaks_custom <- c(-0.2, -0.09, 0.25, 4, max(gene_matrix_all))
colors_custom <- c("white", "#F9E076", "red", "gray")

# Dibujar heatmap
pdf("HeatmapIntrogresiones_porCepa_CladosNeotropical_yalp-2x.pdf",5,10)
pheatmap(
  gene_matrix_all,
  color = colors_custom,
  breaks = breaks_custom,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  fontsize_row = 7,
  annotation_col = annotation_col,
  annotation_legend = TRUE,
  annotation_names_col = TRUE
)
dev.off()
