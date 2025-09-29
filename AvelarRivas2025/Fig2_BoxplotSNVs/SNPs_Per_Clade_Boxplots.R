#### Boxplots No de SNPS ####
library(ggplot2)
library(tidyr)
library(dplyr)
archivo="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/20240920_Submission/AllStrains_WithSNP_Numbers.csv"
df=read.csv(archivo)

head(df)
dim(df)
df <- df %>% mutate(Phylogenetic_Clade = reorder(Phylogenetic_Clade, NumberOfSNPsMDS, FUN = median))

df <- df %>%  group_by(Phylogenetic_Clade) %>% mutate(n_samples = n()) %>% ungroup()

# Step 2: Create a new label for Phylogenetic_Clade including (n=X)
df <- df %>%
  mutate(Phylogenetic_Clade_with_n = paste0(Phylogenetic_Clade, " (n=", n_samples, ")")) %>%
  mutate(Phylogenetic_Clade_with_n = reorder(Phylogenetic_Clade_with_n, NumberOfSNPsMDS, FUN = median))


# Step 3: Plot with reordered labels and count displayed
pdf("250521_NumberOfSNPsMDS.pdf")
ggplot(df, aes(x = Phylogenetic_Clade_with_n, y = NumberOfSNPsMDS)) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Phylogenetic Clade", y = "Number of SNPs MDS") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8 ))   # Rotate x-axis labels
dev.off()

pdf("NumberOfSNPS_SACE.pdf")
ggplot(df, aes(x = Phylogenetic_Clade_with_n, y = NumberOfSNPS_SACE)) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Phylogenetic Clade", y = "NumberOfSNPS_SACE") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8 ))   # Rotate x-axis labels
dev.off()

pdf("NumberOfSNPS_CONC.pdf")
ggplot(df, aes(x = Phylogenetic_Clade_with_n, y = NumberOfSNPS_CONC)) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Phylogenetic Clade", y = "NumberOfSNPS_CONC") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8 ))   # Rotate x-axis labels
dev.off()



ordered_clades <- levels(df$Phylogenetic_Clade_with_n)
indices <- c(36, 1:26, 29, 33:35, 27, 28, 30, 31, 32)
print(ordered_clades[indices])


df$Phylogenetic_Clade_with_n <- factor(df$Phylogenetic_Clade_with_n, 
                                       levels = df$Phylogenetic_Clade_with_n[indices])

ggplot(df, aes(x = Phylogenetic_Clade_with_n, y = NumberOfSNPsMDS)) +
  geom_boxplot() +  # Boxplot with color fill
  geom_point()+
  labs(x = "Phylogenetic Clade", y = "Number of SNPs MDS") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))  # Rotate x-axis labels



clados_primero <- c( "SAM2 (n=4)",  "10. French Guiana human (n=31)",  "Mexican_Agave_1 (n=20)",
    "WB3 (n=3)",  "Tequila_Distillery (n=15)","Mexican_Agave_2 (n=179)")

# 2. Calcular medianas de los otros clados
otros_clados_ordenados <- df %>%
  filter(!Phylogenetic_Clade_with_n %in% clados_primero) %>%
  mutate(Phylogenetic_Clade_with_n = as.character(Phylogenetic_Clade_with_n)) %>%  # <- esta línea es clave
  group_by(Phylogenetic_Clade_with_n) %>%
  summarize(mediana = median(NumberOfSNPsMDS, na.rm = TRUE)) %>%
  arrange(mediana) %>%
  pull(Phylogenetic_Clade_with_n)

# 3. Crear orden final de niveles
orden_final <- c(clados_primero, otros_clados_ordenados)

# 4. Reordenar la variable como factor
df$Phylogenetic_Clade_with_n <- factor(df$Phylogenetic_Clade_with_n, levels = orden_final)

# 5. Gráfico con boxplot + jitter
install.packages("forcats")
#library(forcats)
df <- df %>%
  mutate(Phylogenetic_Clade_with_n = fct_rev(Phylogenetic_Clade_with_n))

#### Figura artículo Revision ####
pdf("Boxplot_NumberSNPs_MDS_data.pdf",5,10)
ggplot(df, aes(x = Phylogenetic_Clade_with_n, y = NumberOfSNPsMDS)) +
  geom_boxplot(outlier.shape = NA) +  # Oculta outliers para no duplicar con jitter
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, color = "red", size = 0.4) +  # línea roja de la mediana
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +  # Puntos dispersos
  labs(x = "Phylogenetic Clade", y = "Number of SNPs MDS") +
  theme_minimal() +
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dev.off()


# hay que:
# Acomodar bien nombres 
# 1-las de french-dairy/alpechin
# 2-NorthAmerican Oak
# 3-SAM1?




#######
archivo="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/20240920_Submission/AllStrains_WithSNP_Numbers.csv"
full_sup_tab=read.csv(archivo)
df=read.csv(archivo)
dim(df)
head(df)
df <- df %>% mutate(Phylogenetic_Clade = reorder(Phylogenetic_Clade, NumberOfSNPsMDS, FUN = median))
df <- df %>%  group_by(Phylogenetic_Clade) %>% mutate(n_samples = n()) %>% ungroup()
# Step 2: Create a new label for Phylogenetic_Clade including (n=X)
df <- df %>%
  mutate(Phylogenetic_Clade_with_n = paste0(Phylogenetic_Clade, " (n=", n_samples, ")")) %>%
  mutate(Phylogenetic_Clade_with_n = reorder(Phylogenetic_Clade_with_n, NumberOfSNPsMDS, FUN = median))

unique_clades <- unique(df$Phylogenetic_Clade_with_n)
indices=c(8,7,3,2,1,5,12,23,21,18,24,10,11,33,17,4,31,20,16,27,19,9,36,32,6,29,25,13,15,14,22,35,28,34,30,26)
indices=c(8,7,3,10,2,1,5,13,24,22,19,25,11,12,34,18,4,32,21,17,28,20,9,37,33,6,30,26,14,16,15,23,36,29,35,31,27)

ordered_levels <- unique_clades[indices]
df$Phylogenetic_Clade_with_n <- factor(df$Phylogenetic_Clade_with_n, levels = ordered_levels)

pdf("NumberOfSNPsMDS_3.pdf",8,4)
ggplot(df, aes(x = Phylogenetic_Clade_with_n, y = NumberOfSNPsMDS)) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Phylogenetic Clade", y = "Number of SNPs MDS") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))  # Rotate x-axis labels
dev.off()


library(dplyr)
library(ggplot2)

# Calcular varianza e IQR por clado
summary_df <- df %>%
  group_by(Phylogenetic_Clade_with_n) %>%
  summarise(
    Variance = var(NumberOfSNPsMDS, na.rm = TRUE),
    IQR = IQR(NumberOfSNPsMDS, na.rm = TRUE)
  )

ggplot(summary_df, aes(x = Phylogenetic_Clade_with_n, y = Variance)) +
  geom_col(fill = "steelblue") +
  labs(x = "Phylogenetic Clade", y = "Variance of Number of SNPs (MDS)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

ggplot(summary_df, aes(x = Phylogenetic_Clade_with_n, y = IQR)) +
  geom_col(fill = "darkorange") +
  labs(x = "Phylogenetic Clade", y = "IQR of Number of SNPs (MDS)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))



pdf("NumberOfSNPS_SACE2.pdf",8,4)
ggplot(df, aes(x = Phylogenetic_Clade_with_n, y = NumberOfSNPS_SACE)) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Phylogenetic Clade", y = "NumberOfSNPS_SACE") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))  # Rotate x-axis labels
dev.off()

pdf("NumberOfSNPS_CONC2.pdf",8,4)
ggplot(df, aes(x = Phylogenetic_Clade_with_n, y = NumberOfSNPS_CONC)) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Phylogenetic Clade", y = "NumberOfSNPS_CONC") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))  # Rotate x-axis labels
dev.off()


#########
n_samples <- df %>%
  filter(!is.na(NumberOfSNPsMDS)) %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(n_samples = n(), .groups = 'drop')

df <- df %>%
  left_join(n_samples, by = "Phylogenetic_Clade") %>%
  mutate(Phylogenetic_Clade_with_n = paste0(Phylogenetic_Clade, " (n=", n_samples, ")")) %>%
  mutate(Phylogenetic_Clade_with_n = reorder(Phylogenetic_Clade_with_n, NumberOfSNPsMDS, FUN = median))
####
# empezar de cero. SOLO MDS####

archivo="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/20240920_Submission/AllStrains_WithSNP_Numbers.csv"
df=read.csv(archivo)

# Contar muestras no NA para cada categoría
n_samples <- df %>%
  filter(!is.na(NumberOfSNPsMDS)) %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(n_samples = n(), .groups = 'drop')
# Calcular la mediana de NumberOfSNPS_SACE para cada categoría
medians <- df %>%
  filter(!is.na(NumberOfSNPsMDS)) %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(median_SACE = median(NumberOfSNPsMDS, na.rm = TRUE), .groups = 'drop')
# Reordenar las categorías manteniendo las categorías específicas al principio
specific_categories <- c("SAM2", "10. French Guiana human", "Mexican_Agave_1", "Tequila_Distillery", "Mexican_Agave_2")
remaining_categories <- setdiff(medians$Phylogenetic_Clade, specific_categories)
ordered_categories <- c(specific_categories, remaining_categories[order(medians$median_SACE[match(remaining_categories, medians$Phylogenetic_Clade)])])
# Asegúrate de que la columna Phylogenetic_Clade esté en el orden correcto
df$Phylogenetic_Clade <- factor(df$Phylogenetic_Clade, levels = ordered_categories)
# Crear el boxplot
pdf("NumberOfSNPsMDS4.pdf",8,4)
ggplot(df, aes(x = Phylogenetic_Clade, y = NumberOfSNPsMDS)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.1) +  # Mostrar todos los puntos
  geom_text(data = n_samples, 
            aes(x = Phylogenetic_Clade, 
                y = max(df$NumberOfSNPsMDS, na.rm = TRUE) + 0.5, 
                label = paste0("n=", n_samples)), 
            size = 3, angle = 90, 
            position = position_nudge(y = 20)) + # Ajustar la posición del texto
  labs(x = "Phylogenetic Clade", y = "NumberOfSNPsMDS") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dev.off()

# empezar de cero. SOLO mapped_to_SACE####

archivo="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/20240920_Submission/AllStrains_WithSNP_Numbers.csv"
df=read.csv(archivo)

# Contar muestras no NA para cada categoría
n_samples <- df %>%
  filter(!is.na(NumberOfSNPS_SACE)) %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(n_samples = n(), .groups = 'drop')
# Calcular la mediana de NumberOfSNPS_SACE para cada categoría
medians <- df %>%
  filter(!is.na(NumberOfSNPS_SACE)) %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(median_SACE = median(NumberOfSNPS_SACE, na.rm = TRUE), .groups = 'drop')
# Reordenar las categorías manteniendo las categorías específicas al principio
specific_categories <- c("SAM2", "10. French Guiana human", "Mexican_Agave_1", "WB3","Mexican_Agave_2", "Tequila_Distillery")
remaining_categories <- setdiff(medians$Phylogenetic_Clade, specific_categories)
#ordered_categories <- c(specific_categories, remaining_categories[order(medians$median_SACE[match(remaining_categories, medians$Phylogenetic_Clade)])])
# Asegúrate de que la columna Phylogenetic_Clade esté en el orden correcto
df$Phylogenetic_Clade <- factor(df$Phylogenetic_Clade, levels = ordered_categories)
pdf("NumberOfSNPS_SACE3.pdf",8,4)
ggplot(df, aes(x = Phylogenetic_Clade, y = NumberOfSNPS_SACE)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = .52, alpha = 0.1) +  # Mostrar todos los puntos
  geom_text(data = n_samples, 
            aes(x = Phylogenetic_Clade, 
                y = max(df$NumberOfSNPS_SACE, na.rm = TRUE) + 0.5, 
                label = paste0("n = ", n_samples)), 
            size = 3, angle = 90, 
            position = position_nudge(y = 20)) + # Ajustar la posición del texto
  labs(x = "Phylogenetic Clade", y = "Number of SNPs SACE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dev.off()
# empezar de cero. SOLO mapped_to_CONC####

archivo="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/20240920_Submission/AllStrains_WithSNP_Numbers.csv"
df=read.csv(archivo)

# Contar muestras no NA para cada categoría
n_samples <- df %>%
  filter(!is.na(NumberOfSNPS_CONC)) %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(n_samples = n(), .groups = 'drop')
# Calcular la mediana de NumberOfSNPS_SACE para cada categoría
medians <- df %>%
  filter(!is.na(NumberOfSNPS_CONC)) %>%
  group_by(Phylogenetic_Clade) %>%
  summarise(median_SACE = median(NumberOfSNPS_CONC, na.rm = TRUE), .groups = 'drop')
# Reordenar las categorías manteniendo las categorías específicas al principio
specific_categories <- c("SAM2", "10. French Guiana human", "Mexican_Agave_1", "WB3","Mexican_Agave_2", "Tequila_Distillery")
remaining_categories <- setdiff(medians$Phylogenetic_Clade, specific_categories)
ordered_categories <- c(specific_categories, remaining_categories[order(medians$median_SACE[match(remaining_categories, medians$Phylogenetic_Clade)])])
# Asegúrate de que la columna Phylogenetic_Clade esté en el orden correcto
df$Phylogenetic_Clade <- factor(df$Phylogenetic_Clade, levels = ordered_categories)
pdf("NumberOfSNPS_CONC3.pdf",8,4)
ggplot(df, aes(x = Phylogenetic_Clade, y = NumberOfSNPS_CONC)) +
  geom_boxplot() +
  geom_text(data = n_samples, 
            aes(x = Phylogenetic_Clade, 
                y = max(df$NumberOfSNPS_SACE, na.rm = TRUE) + 0.5, 
                label = paste0("n = ", n_samples)), 
            size = 3, angle = 90, 
            position = position_nudge(y = 20)) + # Ajustar la posición del texto
  labs(x = "Phylogenetic Clade", y = "NumberOfSNPS_CONC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dev.off()

# hacer SACE y CONC en el mismo plot #####
library(dplyr)

df <- df %>%
  group_by(Phylogenetic_Clade) %>%
  mutate(
    n_valid_CONC = sum(!is.na(NumberOfSNPS_CONC)),  # cuenta solo valores válidos
    Phylogenetic_Clade_with_n = paste0(Phylogenetic_Clade, " (n=", n_valid_CONC, ")")
  ) %>%
  ungroup() %>%
  mutate(
    Phylogenetic_Clade_with_n = reorder(Phylogenetic_Clade_with_n, NumberOfSNPS_CONC, FUN = median)
  )


clados_primero <- c(   "Mexican_Agave_1 (n=20)",
                     "WB3 (n=3)",  "Tequila_Distillery (n=15)","Mexican_Agave_2 (n=179)","SAM2 (n=4)",  "10. French Guiana human (n=29)")
otros_clados_ordenados <- df %>%
  filter(!Phylogenetic_Clade_with_n %in% clados_primero) %>%
  mutate(Phylogenetic_Clade_with_n = as.character(Phylogenetic_Clade_with_n)) %>%  # <- esta línea es clave
  group_by(Phylogenetic_Clade_with_n) %>%
  summarize(mediana = median(NumberOfSNPS_CONC, na.rm = TRUE)) %>%
  arrange(mediana) %>%
  pull(Phylogenetic_Clade_with_n)
orden_final <- c(clados_primero, otros_clados_ordenados)
df$Phylogenetic_Clade_with_n <- factor(df$Phylogenetic_Clade_with_n, levels = orden_final)


df_long <- df %>%
  pivot_longer(cols = c(NumberOfSNPS_SACE, NumberOfSNPS_CONC), 
               names_to = "Variable", values_to = "Value")

# Contar muestras no NA para cada categoría y cada variable
n_samples <- df_long %>%
  filter(!is.na(Value)) %>%
  group_by(Phylogenetic_Clade, Variable) %>%
  summarise(n_samples = n(), .groups = 'drop')

# Crear el boxplot con las dos variables diferenciadas por colorv
#### FIGURA ARTÍCULO CONC VS SACE
pdf("250608_SACEyCONC.pdf",8,4)
ggplot(df_long, aes(x = Phylogenetic_Clade_with_n, y = Value, fill = Variable)) +
  geom_boxplot(width = 0.75, position = position_dodge(width = 0.8)) +  # Mismo ancho que dodge
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 0.75, alpha = 0.25) +  # puntos con color opcional
  labs(x = "Phylogenetic Clade", y = "Number of SNPs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "top")
dev.off()





## otra iteracion incluyendo y excluyendo introgresiones####
library(dplyr)
library(tidyr)
library(ggplot2)

# Especificar el orden deseado para las categorías
ordered_categories <- c("1. Wine/European", "2. Alpechin", "10. French Guiana human", 
                        "SAM2", "Tequila_Distillery", 
                        "Mexican_Agave_1", "Mexican_Agave_2", "17. Taiwanese")

# Convertir la columna 'Phylogenetic_Clade' en factor con el orden especificado
df_long$Phylogenetic_Clade <- factor(df_long$Phylogenetic_Clade, levels = ordered_categories)

# Contar muestras no NA para cada categoría y cada variable
n_samples <- df_long %>%
  filter(!is.na(Value)) %>%
  group_by(Phylogenetic_Clade, Variable) %>%
  summarise(n_samples = n(), .groups = 'drop')
pdf("IntvsNOINT.pdf")
ggplot(df_long, aes(x = Phylogenetic_Clade, y = Value, fill = Variable)) +
  geom_boxplot(width = 0.75, position = position_dodge(width = 0.8)) +  # Boxplots más delgados con posición separada
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.5, alpha = 0.4) +  # Puntos desplazados ligeramente
  geom_text(data = n_samples, 
            aes(x = Phylogenetic_Clade, 
                y = max(df_long$Value, na.rm = TRUE) + 0.5, 
                label = paste0("n = ", n_samples), fill = Variable), 
            size = 3, angle = 90, 
            position = position_nudge(y = 20)) + # Ajustar la posición del texto
  labs(x = "Phylogenetic Clade", y = "Number of SNPs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "none")  # Quitar la leyenda
dev.off()



clados_base <- c("1. Wine/European","2. Alpechin", "Mexican_Agave_1", "WB3", "Tequila_Distillery",
                 "Mexican_Agave_2", "SAM2", "10. French Guiana human",
                  "17. Taiwanese")
df$Phylogenetic_Clade <- as.character(df$Phylogenetic_Clade)
df$Phylogenetic_Clade_with_n <- as.character(df$Phylogenetic_Clade_with_n)

df_mod <- df %>%
  mutate(Clado_agrupado = ifelse(
    Phylogenetic_Clade %in% clados_base,
    Phylogenetic_Clade_with_n,  # nombre completo con n
    "others"
  ))
orden_clados_plot <- c(
  df %>%
    filter(Phylogenetic_Clade %in% clados_base) %>%
    distinct(Phylogenetic_Clade, Phylogenetic_Clade_with_n) %>%
    arrange(match(Phylogenetic_Clade, clados_base)) %>%
    pull(Phylogenetic_Clade_with_n),
  "others"
)

df_mod$Clado_agrupado <- factor(df_mod$Clado_agrupado, levels = orden_clados_plot)
df_mod_long <- df_mod %>%
  pivot_longer(cols = c(NumberOfSNPS_SACE, NumberOfSNPS_CONC), 
               names_to = "Variable", values_to = "Value")

pdf("IntVsNoInt_onlyfewClades2.pdf",4,4)
ggplot(df_mod_long, aes(x = Clado_agrupado, y = Value, fill = Variable)) +
  geom_boxplot(width = 0.75, position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             size = 0.75, alpha = 0.25) +
  labs(x = "Phylogenetic Clade", y = "Number of SNPs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "top")
dev.off()




##### las diferencias ####

# Crear una nueva columna con las diferencias entre NumberOfSNPS_SACE y NumberOfSNPS_CONC
df <- df %>%
  mutate(Difference = NumberOfSNPS_SACE - NumberOfSNPS_CONC)

# Convertir la columna 'Phylogenetic_Clade' en factor con el orden especificado
df$Phylogenetic_Clade <- factor(df$Phylogenetic_Clade, levels = ordered_categories)

# Crear el boxplot de las diferencias
ggplot(df, aes(x = Phylogenetic_Clade, y = Difference)) +
  geom_boxplot(width = 0.75) +  # Boxplots más delgados
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.75, alpha = 0.25) +  # Mostrar todos los puntos con jitter
  labs(x = "Phylogenetic Clade", y = "Difference (SACE - CONC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))


ggplot(df, aes(x = NumberOfSNPS_SACE, y = NumberOfSNPS_CONC, color = Phylogenetic_Clade)) +
  geom_point(size = 2, alpha = 0.8) +  # Puntos
#  scale_color_manual(values = palette_colors) +  # Colores manuales
  labs(x = "Number of SNPs SACE", y = "Number of SNPs CONC") +
  theme_minimal() +
  theme(#legend.position = "none",  # Sin leyenda
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8))




### Febrero 2025, solo la matriz de Peter et al #####
archivo1="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/figs_code/Snps_perSample_1011Matrix_metadata.csv"
names=read.csv(archivo1)
dim(names)


names <- names %>%
  group_by(Clades) %>%
  mutate(Clades_n = paste0(Clades, " (n=", n(), ")")) %>%
  ungroup()

ggplot(names, aes(x = reorder(Clades_n, Total.number.of.SNPs, FUN = median), y = Total.number.of.SNPs)) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Clades", y = "Number of SNPs Peter_et_al_2018") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8 ))   # Rotate x-axis labels

ggplot(names, aes(x = reorder(Clades_n, Total.number.of.SNPs, FUN = median), y = SNPs )) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Clades", y = "Number of SNPs SNPs from Matrix1011") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8 ))   # Rotate x-axis labels


names[names$Clades=="9.2 other neotropical",]

## las de victor 3034#####
archivo1="C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/figs_code/SNPS_PerSample_Loegler2024.txt .csv"
names=read.csv(archivo1)
dim(names)
head(names)

names <- names %>%
  group_by(Clade) %>%
  mutate(Clades_n = paste0(Clade, " (n=", n(), ")")) %>%
  ungroup()

ggplot(names, aes(x = reorder(Clades_n, SNPs, FUN = median), y = SNPs)) +
  geom_boxplot() +  # Boxplot with color fill
  labs(x = "Clades", y = "Number of SNPs Loegler et al 2024") +  # Labels
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8 ))   # Rotate x-axis labels




#Marzo 2025 usando vcfR
install.packages("vcfR")
library(vcfR)
miarchivo="C:/Users/HP/Dropbox/Posdoc/Supplementary Files From_Peter_et_al_2018/Loegler_2024/GVCF_3039samples.tar.gz"
vcf <- read.vcfR(miarchivo)
gt <- extract.gt(vcf)

snp_counts <- apply(gt, 2, function(x) sum(x %in% c("0/1", "1/1", "1|0", "1|1")))
print(snp_counts)



#### PI  ####
library(ggplot2)
library(dplyr)
library(readr)

# Define el directorio donde están los archivos .pi
dir_path <- "c:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/figs_code/PopParams_SACE1264/Clades"

# Lista de archivos que terminan en .pi
pi_files <- list.files(dir_path, pattern = "\\.pi$", full.names = TRUE)

# Leer todos los archivos válidos
pi_data_list <- lapply(pi_files, function(file) {
  # Intenta leer el archivo, atrapando errores
  df <- tryCatch(read_table(file, col_names = TRUE), error = function(e) NULL)
  
  # Verifica que se haya leído bien y que tenga la columna "PI"
  if (!is.null(df) && "PI" %in% names(df)) {
    return(data.frame(PI = df$PI, File = basename(file)))
  } else {
    warning(paste("Archivo inválido o sin columna PI:", basename(file)))
    return(NULL)
  }
})

# Eliminar elementos NULL del resultado
pi_data_list <- pi_data_list[!sapply(pi_data_list, is.null)]

# Verifica si hay datos antes de graficar
if (length(pi_data_list) == 0) {
  stop("Ningún archivo válido con columna 'PI' fue encontrado.")
}

# Combinar todos los data.frames
pi_data <- bind_rows(pi_data_list)

pi_data$File <- gsub("_Window_results.windowed.pi$", "", pi_data$File)

# Calcular media y mediana por archivo
summary_stats <- pi_data %>%
  group_by(File) %>%
  summarise(Media = mean(PI, na.rm = TRUE),
            Mediana = median(PI, na.rm = TRUE)) %>%
  tidyr::pivot_longer(cols = c(Media, Mediana), names_to = "Estadística", values_to = "PI")

# Crear el gráfico con solo media y mediana
ggplot(summary_stats, aes(x = File, y = PI, fill = Estadística)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +
  labs(
       x = "Clade",
       y = "Pi") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_blank())

#acomodar de menos a mas
#poner a las neotropical al principio
#agregar n
