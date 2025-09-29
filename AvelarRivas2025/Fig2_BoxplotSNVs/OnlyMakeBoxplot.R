

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