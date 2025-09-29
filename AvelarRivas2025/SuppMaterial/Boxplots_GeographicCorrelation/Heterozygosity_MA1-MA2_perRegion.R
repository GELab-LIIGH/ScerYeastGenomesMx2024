# BOXPLOT HETEROZYGOSITY PER SUBPOPLATION PER STATE
archivo="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/20240920_Submission/AllStrains_WithSNP_Numbers.csv"
full_sup_tab=read.csv(archivo)
colnames(full_sup_tab)
archivo="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/20241111_NatComm_Submission/AgaveYeast_SupplementaryTables_MOD.csv"
mx_metadata=read.csv(archivo)
colnames(mx_metadata)


# Agregar Region desde mx_metadata
full_sup_tab <- merge(full_sup_tab,
                      mx_metadata[, c("Sample_Name", "Region")],
                      by = "Sample_Name",
                      all.x = TRUE)

# Filtrar solo los que están en Mexican_Agave_1 o _2
subset_data <- subset(full_sup_tab, 
                      Phylogenetic_Clade %in% c("Mexican_Agave_1", "Mexican_Agave_2") & 
                        !is.na(Region) & !is.na(Heterozygosity))

# Asegurar que Heterozygosity sea numérica
subset_data$Heterozygosity <- as.numeric(as.character(subset_data$Heterozygosity))

# Hacer el boxplot
boxplot(Heterozygosity ~ Region, data = subset_data,
        main = "Heterozygosity por Región",
        xlab = "Región", ylab = "Heterozygosidad",
        col = rainbow(length(unique(subset_data$Region))),
        las = 2, cex.axis = 0.8)

ggplot(subset_data, aes(x = Region, y = Heterozygosity, fill = Region)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Heterozygosity por Región",
       x = "Región", y = "Heterozygosidad") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


region_order <- c("Northwest", "West-I", "West-II", 
                  "CentralHighlands", "BalsasBasin", 
                  "SouthCentral", "Northeast")
region_order <- c("Northeast", "SouthCentral", "BalsasBasin", 
                  "CentralHighlands", "West-II", "West-I", "Northwest")
# Filtrado como antes
subset_data <- subset(full_sup_tab, 
                      Phylogenetic_Clade %in% c("Mexican_Agave_1", "Mexican_Agave_2") &
                        !is.na(Region) & !is.na(Heterozygosity))

# Convertir Region a factor con el orden deseado
subset_data$Region <- factor(subset_data$Region, levels = region_order)

# Graficar
# 1. Cepas Mexican_Agave_1
ma1 <- subset(full_sup_tab, Phylogenetic_Clade == "Mexican_Agave_1")

# 2. Cepas que no están en mx_metadata
missing_samples <- setdiff(ma1$Sample_Name, mx_metadata$Sample_Name)

# 3. Subconjunto de full_sup_tab con esas cepas faltantes
missing_data <- subset(ma1, Sample_Name %in% missing_samples)

# 4. Asignar "Northeast" como región
missing_data$Region <- "Northeast"

# 5. Agregar esos datos al subset_data existente
subset_data_extended <- rbind(subset_data, missing_data)



ggplot(subset_data, aes(x = Region, y = Heterozygosity, fill = Region)) +
  #geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.35) +
  coord_flip() +
  theme_minimal() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(title = "Heterozygosity per Region (only MA1 and MA2)",
       x = "Región", y = "Heterozygosidad") +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))

pdf("Heterozygosity_MA1-MA2.pdf",5,5)
ggplot(subset_data_extended, aes(x = Region, y = Heterozygosity, color = Region)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Only MA1 and MA2",
       x = "Región", y = "Heterozygosity") +
  scale_color_brewer(palette = "Set3") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
dev.off()
