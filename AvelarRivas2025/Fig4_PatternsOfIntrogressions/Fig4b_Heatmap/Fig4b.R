#Viene de figs_code/figsR/AnalyzeIntrogressionTrees.R L1144
library(readr)
merged_final_combined_table_2 <- read_csv("final_combined_table_2_perGene_250609.csv")


merged_genes_expanded <- merged_final_combined_table_2 %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = trimws(Genes))  # Elimina espacios extra si los hay
dim(merged_genes_expanded)
View(merged_genes_expanded)

merged_genes_expanded$Gene <- substr(merged_genes_expanded$Genes, 
                                     nchar(merged_genes_expanded$Genes) - 7, 
                                     nchar(merged_genes_expanded$Genes))


# hacer heatmap para figura 4b:
library(dplyr)
library(tidyr)

gene_matrix_Neotropical <- merged_genes_expanded %>%
  filter(Phylogenetic_Clade %in% c("SAM2", "WB3","Mexican_Agave_1","Mexican_Agave_2","French Guiana")) %>%  # Cepas de SAM2 o WB3
  select(Strain, Genes) %>%                 # Tomar solo columnas relevantes
  filter(!is.na(Genes)) %>%                 # Remover genes NA
  distinct() %>%                            # Evitar duplicados Strain-Gene
  mutate(presence = 1) %>%                  # Indicador de presencia
  pivot_wider(names_from = Genes, values_from = presence, values_fill = 0) %>%  # Expande a matriz
  arrange(Strain)

library(pheatmap)
gene_matrix_for_heatmap <- gene_matrix_Neotropical %>%
  tibble::column_to_rownames("Strain")


gene_matrix_clade <- merged_genes_expanded %>%
  filter(Phylogenetic_Clade %in% c("SAM2", "WB3", "Tequila_Distillery", "Mexican_Agave_1", "Mexican_Agave_2", "French Guiana","Alpechin")) %>%
  select(Strain, Phylogenetic_Clade, Genes) %>%
  filter(!is.na(Genes)) %>%
  distinct() %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Genes, values_from = presence, values_fill = 0)

# Paso 2: Calcular porcentaje de presencia por clado
gene_percent_by_clade <- gene_matrix_clade %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(across(where(is.numeric), ~mean(.x) * 100)) %>%
  as.data.frame()

rownames(gene_percent_by_clade) <- gene_percent_by_clade$Phylogenetic_Clade
gene_percent_by_clade$Phylogenetic_Clade <- NULL

color_palette <- colorRampPalette(c("grey90", "blue"))(100)


orden_filas <- c("Mexican_Agave_1", "Mexican_Agave_2", "Tequila_Distillery", 
                 "WB3", "SAM2", "French Guiana", "Alpechin")

# Reordenar las filas manualmente
gene_percent_by_clade <- gene_percent_by_clade[orden_filas, ]

# Definir los cortes y colores
breaks <- c(0, 0.001, seq(1, 100, length.out = 99))  # 0 exclusivo → blanco
colors <- c("white", colorRampPalette(c("grey80", "blue"))(99))

# Dibujar el heatmap
pdf("Heatmap_Fig3b_PresenciaAusenciaGenes2.pdf", 6, 7)
pheatmap(
  gene_percent_by_clade,
  color = colors,
  breaks = breaks,
  cluster_rows = FALSE,  # No reordenar filas automáticamente
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10
)
dev.off()