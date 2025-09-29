library(ape)
library(dplyr)
library(stringr)
  
get_descendants <- function(tree, node) {
  descendants <- c()
  nodes_to_check <- node  # Inicializar con el nodo actual
  
  while (length(nodes_to_check) > 0) {
    new_nodes <- tree$edge[tree$edge[,1] %in% nodes_to_check, 2]  # Hijos directos
    descendants <- c(descendants, new_nodes)
    nodes_to_check <- new_nodes[new_nodes > length(tree$tip.label)]  # Solo nodos internos
  }
  
  # Retornar solo los tips
  return(descendants[descendants <= length(tree$tip.label)])
}
extract_introgression <- function(filename) {
  matches <- regmatches(filename, regexec("Matrix_YPS138_([0-9A-Z]+)_//d+_YPS138_([0-9A-Z]+)_filtered", filename))
  if (length(matches[[1]]) == 3) {
    return(list(start = matches[[1]][2], end = matches[[1]][3]))
  } else {
    return(list(start = NA, end = NA))  # Manejo de errores si no coincide
  }
}
analyze_tree <- function(treefile, reference_strain) {
  tree <- read.tree(treefile)
  #print(treefile)
  #tree <- read.tree(treefiles[9])
  name_replacements <- c(
    "JS445c1SAPA" = "WINE-JS445c1SAPA",
    "JS611c1SAPA" = "Taiwan-JS611c1SAPA",
    "SRR4074385" = "SpA-SRR4074385",
    "SRR4074411" = "SpB-YPS138-SRR4074411",
    "SRR4074412" = "SpB-Bra-SRR4074412",
    "YMX005537" = "SpB-Tank-YMX005537",
    "SRR7500262" = "SpC-SRR7500262"
  )
  
  tree$tip.label <- ifelse(tree$tip.label %in% names(name_replacements),
                           name_replacements[tree$tip.label],
                           tree$tip.label)
  
  #plot(tree)
  #plot(tree, type = "unrooted", use.edge.length = TRUE, show.tip.label = TRUE, no.margin = TRUE, cex = 0.75, main = file)
  #nodelabels()
  #best_node <- find_best_sp_node(tree, reference_strain)
  #tree <- root(tree, node = best_node, resolve.root = TRUE)
  
  # Re-root the tree at Taiwan-JS611c1SAPA
  if ("Taiwan-JS611c1SAPA" %in% tree$tip.label) {
    tree <- root(tree, outgroup = "Taiwan-JS611c1SAPA", resolve.root = TRUE)
  }
  # Find the node closest to the reference strain
  ref_node <- which(tree$tip.label == reference_strain)
  parent_node <- tree$edge[which(tree$edge[,2] == ref_node), 1]
  # Compute distances to all nodes
  dist_nodes <- cophenetic(tree)[reference_strain, ]
  # Get all descendants of the closest node
  all_child_nodes <- get_descendants(tree, parent_node)
  child_labels <- setdiff(tree$tip.label[all_child_nodes], reference_strain)  # Remove reference strain
  grandparent_node <- tree$edge[which(tree$edge[,2] == parent_node), 1]
  if (length(grandparent_node) == 0) {
    grandparent_node <- NA  # Manejar el caso donde el nodo padre es la raíz
  }
  all_grandchildren <- get_descendants(tree, grandparent_node)
  extra_labels <- setdiff(tree$tip.label[all_grandchildren], c(child_labels, reference_strain))  # Remove ref strain from extra
  
  # Extract introgression start and end
  introgression_info <- extract_introgression(basename(treefile))
  
  # Convert distance matrix to named vector (each node as a column)
  dist_data <- as.data.frame(t(dist_nodes))

  sorted_distances <- sort(dist_nodes[names(dist_nodes) != reference_strain])  # Exclude self
  closest_relative <- names(sorted_distances)[1]
  second_closest_relative <- names(sorted_distances)[2]
  
  if (closest_relative == ifelse(length(child_labels) > 0, paste(child_labels, collapse = "."), NA) ){
    IntrOrigin_strict = closest_relative
  }
  else{
      IntrOrigin_strict = NA
      }
  
  expected_strains <- name_replacements  # Este es un named vector: reemplazos
  for (strain in expected_strains) {
    if (!(strain %in% colnames(dist_data))) {
      dist_data[[strain]] <- NA  # Añadir columna con NA si falta
    }
  }
  
  # Reordenar columnas en orden estándar
  dist_data <- dist_data[, sort(names(dist_data))]
  
  
  # Create final output as a single row
  result <- data.frame(
    Intr_Start = introgression_info$start,
    Intr_End = introgression_info$end,
    Strain = reference_strain,
    LowestDistance = closest_relative,
    SecondLowestDistance = second_closest_relative,
    Phylogeny_Siblings = ifelse(length(child_labels) > 0, paste(child_labels, collapse = "."), NA),
    Itr_Origin = IntrOrigin_strict,
    Cousins = ifelse(length(extra_labels) > 0, paste(extra_labels, collapse = "."), NA),
    File = treefile,
    dist_data  
  )
  #result <- cbind(result, dist_data)  # Merge with distance data
  return(result)
}
is_true_introgression <- function(siblings, cousins, lowest_distance, second_lowest_distance) {
  # thisrow=10
  # siblings=results_df$Phylogeny_Siblings[thisrow]
  # cousins=results_df$Cousins[thisrow]
  # lowest_distance=results_df$LowestDistance[thisrow]
  # second_lowest_distance=results_df$SecondLowestDistance[thisrow]
   
  combined <- paste(siblings, cousins, sep = ".")  # Unir las dos columnas
  sapa_count <- str_count(combined, "SAPA")        # Contar ocurrencias de "SAPA"
  
  # Verificar si LowestDistance o SecondLowestDistance contienen "SAPA"
  has_sapa_in_distances <- grepl("SAPA", lowest_distance) & grepl("SAPA", second_lowest_distance)
  
  # True si hay menos de 2 ocurrencias de "SAPA" y no está en las distancias más bajas
  return(sapa_count < 2 & !has_sapa_in_distances)
}
determine_true_origin <- function(true_introgression, lowest_distance, siblings, cousins) {
  if (!true_introgression) return(NA) # Solo evaluamos si TrueIntrogressions es TRUE
  
  # Verificar si LowestDistance y Siblings son exactamente iguales
  if (lowest_distance == siblings) return(lowest_distance)
  
  combined <- unique(unlist(str_split(paste(siblings, cousins, sep = "."), "//."))) # Separar valores
  
  # Casos específicos corregidos
  if (all(sapply(combined, function(x) str_starts(x, "SpB")))) return("SpB_Undefined")
  if (all(sapply(combined, function(x) str_starts(x, "SpB") | str_starts(x, "SpC")))) return("Sp_American")
  if (all(sapply(combined, function(x) str_starts(x, "SpA") | str_starts(x, "Wine")))) return("SpA")
  
  return("Sp_undefined")
}
find_best_sp_node <- function(tree, reference_strain) {
  # Function to get descendant tips of a node
  get_descendants <- function(tree, node) {
    return(tree$tip.label[phytools::getDescendants(tree, node)])
  }
  
  best_node <- NULL
  max_sp_count <- -1
  min_sapa_count <- Inf
  
  for (node in (length(tree$tip.label) + 1):max(tree$edge)) {
    descendants <- get_descendants(tree, node)
    
    # Exclude reference strain from counts
    descendants <- setdiff(descendants, reference_strain)
    
    sp_count <- sum(grepl("^Sp", descendants))      # Count Sp* tips
    sapa_count <- sum(grepl("SAPA", descendants))  # Count SAPA tips
    
    # Update best node if conditions are met
    if (sp_count > max_sp_count || (sp_count == max_sp_count && sapa_count < min_sapa_count)) {
      best_node <- node
      max_sp_count <- sp_count
      min_sapa_count <- sapa_count
    }
  }
  
  return(best_node)
}
best_node <- find_best_sp_node(tree, reference_strain)
update_true_introgressions <- function(results_df) {
  for (i in seq_len(nrow(results_df))) {
    
    if( is.na(results_df$TrueOrigin[i]) ){
      next
    }
    # Modificar TrueOrigin si es "Sp_undefined"
    if (results_df$TrueOrigin[i] == "Sp_undefined") {
      siblings <- unlist(str_split(results_df$Phylogeny_Siblings[i], "//."))
      
      if (all(str_starts(siblings, "Sp"))) {
        results_df$TrueOrigin[i] <- "Sp_und_outgroup"
      }
    }
    
    # Si TrueOrigin es "WINE-JS445c1SAPA", modificar usando el árbol
    if (results_df$TrueOrigin[i] == "WINE-JS445c1SAPA") {
      tree_file <- results_df$File[i]
      
      if (!is.na(tree_file) && file.exists(tree_file)) {
        tree <- read.tree(tree_file)
        reference_strain <- results_df$Strain[i]
        
        best_node <- find_best_sp_node(tree, reference_strain)
        
        # Enraizar con el best_node
        rooted_tree <- root(tree, node = best_node, resolve.root = TRUE)
        
        # Recalcular Siblings y Cousins
        ref_node <- which(rooted_tree$tip.label == reference_strain)
        parent_node <- rooted_tree$edge[which(rooted_tree$edge[,2] == ref_node), 1]
        
        get_descendants <- function(tree, node) {
          return(tree$tip.label[phytools::getDescendants(tree, node)])
        }
        
        siblings <- setdiff(get_descendants(rooted_tree, parent_node), reference_strain)
        
        grandparent_node <- rooted_tree$edge[which(rooted_tree$edge[,2] == parent_node), 1]
        cousins <- setdiff(get_descendants(rooted_tree, grandparent_node), c(siblings, reference_strain))
        
        # Unir Siblings y Cousins
        combined <- c(siblings, cousins)
        
        # Si ambas cepas con SAPA están presentes, cambiar TrueIntrogressions a FALSE y TrueOrigin a NA
        if (sum(grepl("SAPA", combined)) >= 2) {
          results_df$TrueIntrogressions[i] <- FALSE
          results_df$TrueOrigin[i] <- NA
        }
      }
    }
  }
  
  return(results_df)
}

base_directory <- "C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/figs_code/SACE_ToNatC/IntrogrIndividualPhylogenes/250409_allPhylogenies"
reference_strains <- list.dirs(base_directory, full.names = TRUE, recursive = FALSE)
analyze_all_strains <- function(reference_strain) {
  print(reference_strain)
  directory <- file.path(base_directory, reference_strain)
  
  # Comprobar si el directorio tiene archivos .treefile
  treefiles <- list.files(directory, pattern = "//.treefile$", full.names = TRUE)
  
  if (length(treefiles) > 0) {
    # Aplicar la función de análisis de árboles a todos los archivos .treefile
    results <- lapply(treefiles, analyze_tree, reference_strain = reference_strain)
    
    # Filtrar los resultados NULL
    results <- Filter(Negate(is.null), results)
    
    if (length(results) > 0) {
      # Convertir los resultados en un solo dataframe
      results_df <- do.call(rbind, results)
      
      # Imprimir el nombre de la cepa de referencia y la tabla de la séptima columna
      print(paste("Referencia:", reference_strain))
      print(table(results_df[, 7]))
      
      # Guardar los resultados en un archivo CSV
      write.csv(results_df, file = paste0(reference_strain, "_results.csv"), row.names = FALSE)
    } else {
      print(paste("No se encontraron resultados para", reference_strain))
    }
  } else {
    print(paste("No se encontraron archivos .treefile en el directorio de", reference_strain))
  }
}
resultsdir=paste(base_directory, "/ResultadosPorCepa/",sep="")
for (reference_strain in reference_strains) {
  strain_name <- basename(reference_strain)
  if (strain_name == "ResultadosPorCepa" || strain_name == "results" ){
    next
  }
  #reference_strain = reference_strains[19]
  strain_name <- basename(reference_strain)
  print(strain_name)
  
  treefiles <- list.files(reference_strain, pattern = "//.treefile$", full.names = TRUE)
  results <- lapply(treefiles, analyze_tree, reference_strain = strain_name)
  results <- Filter(Negate(is.null), results)
  results_df <- do.call(rbind, results)
  
  results_df <- results_df %>%
    mutate(
      TrueIntrogressions = mapply(is_true_introgression, Phylogeny_Siblings, Cousins, LowestDistance, SecondLowestDistance),
      TrueOrigin = mapply(determine_true_origin, TrueIntrogressions, LowestDistance, Phylogeny_Siblings, Cousins)
    )
  
  column_names <- colnames(results_df)
  last_two <- tail(column_names, 2)
  new_order <- c(column_names[1:3], last_two, column_names[4:(length(column_names)-2)])
  results_df <- results_df[, new_order]
  results_df <- update_true_introgressions(results_df)
  
  write.csv(results_df, file = paste0(resultsdir, strain_name, "_results_250409.csv"), row.names = FALSE)
  
  print(table(results_df[,5]))
  print(table(results_df[,9]))
}


# Definir las cepas de interés
strains_of_interest <- c("WINE-JS445c1SAPA", "Taiwan-JS611c1SAPA", "SpA-SRR4074385", 
                         "SpB-YPS138-SRR4074411", "SpB-Bra-SRR4074412", 
                         "SpB-Tank-YMX005537", "SpC-SRR7500262","Sp_American","SpB_Undefined","Sp_undefined", "Sp_und_outgroup")
strains_of_interest=toupper(strains_of_interest)
results_directory <- resultsdir # "C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/figs_code/figsR/"
results_files <- list.files(results_directory, pattern = "_results_250409.csv$", full.names = TRUE)
final_df <- data.frame(matrix(ncol = length(strains_of_interest), nrow = length(results_files)))
colnames(final_df) <- strains_of_interest
rownames(final_df) <- basename(results_files)

# Función para contar las cepas en la columna 7 de cada archivo
count_strains_in_file <- function(file) {
  # Leer el archivo CSV
  results_df <- read.csv(file)
  
  # Asegurarse de que la columna 7 sea de tipo character
  results_df[, 5] <- as.character(results_df[, 5])
  
  # Limpiar los valores de la columna 7 (eliminar espacios adicionales)
  results_df[, 5] <- trimws(results_df[, 5])
  
  # Convertir la columna 7 y las cepas a minúsculas para evitar problemas de mayúsculas/minúsculas
  results_df[, 5] <- toupper(results_df[, 5])
  
  # Usar table() para contar las ocurrencias de cada valor en la columna 7
  strain_table <- table(results_df[, 5])
  
  # Inicializar un vector de conteos para las cepas de interés
  strain_counts <- rep(0, length(strains_of_interest))
  names(strain_counts) <- strains_of_interest
  
  # Asignar los conteos correspondientes de la tabla
  for (strain in strains_of_interest) {
    # Si la cepa está presente en la tabla, asignamos el conteo
    if (strain %in% names(strain_table)) {
      strain_counts[strain] <- strain_table[strain]
    }
  }
  
  return(strain_counts)
}


# Iterar sobre todos los archivos y contar las cepas
for (i in 1:length(results_files)) {
  # Llamar a la función para contar las cepas en cada archivo
  strain_counts <- count_strains_in_file(results_files[i])
  
  # Almacenar los conteos en el dataframe final
  final_df[i, ] <- strain_counts
}
final_df
# Guardar el dataframe final en un archivo CSV
rownames(final_df) <- sub("_.*", "", rownames(final_df))
write.csv(final_df, file = "250409_strain_counts_summary.csv", row.names = TRUE)
head(final_df)
print(final_df)


data <- read.csv("C:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/20241111_NatComm_Submission/SampleSheets/SampleSheet_To_Pick_Randomly.csv")
library(ggplot2)
library(dplyr)
library(tidyr)

# Convertir rownames a columna en final_df
final_df <- final_df %>%
  tibble::rownames_to_column(var = "ID")

# Unir final_df con data para obtener Phylogenetic_Clade
merged_df <- final_df %>%
  left_join(data, by = "ID")



clades_to_plot <- c("Mexican_Agave_1", "Mexican_Agave_2", "Tequila_Distillery","French Guiana",
                    "WB3","SAM2","Alpechin")

# Variables a excluir
excluded_variables <- c("SPC-SRR7500262", "TAIWAN-JS611C1SAPA", "WINE-JS445C1SAPA","SP_UNDEFINED")

ordered_vars <- c(  "SPB-TANK-YMX005537",  "SPB-YPS138-SRR4074411",  "SPB-BRA-SRR4074412",
  "SPB_UNDEFINED",  "SPC-SRR7500262",  "SP_AMERICAN",  "SPA-SRR4074385",  "SP_UND_OUTGROUP")

# Aplicar filtro y convertir a factor ordenado
filtered_data <- summary_stats %>%
  filter(Phylogenetic_Clade %in% clades_to_plot,
         Variable %in% ordered_vars) %>%
  mutate(Variable = factor(Variable, levels = ordered_vars))










































## ####
library(dplyr)
install.packages("readr")
library(readr)

all_blocks <- read_csv("C:/Users/HP/Desktop/FromIvan/All_blocks.csv")

# Inicializar lista para almacenar los resultados combinados
merged_all_results <- list()

# Iterar sobre cada archivo de resultados
for (file in results_files) {
  # Leer archivo de resultados
  result_df <- read_csv(file)
  
  # Crear columna Introg_name
  result_df <- result_df %>%
    mutate(Introg_name = paste0("YPS138_", Intr_Start, "_2_YPS138_", Intr_End))
  
  # Hacer merge con all_blocks por columna "Name"
  merged_df <- result_df %>%
    left_join(all_blocks, by = c("Introg_name" = "Name"))
  
  # Agregar al acumulador
  merged_all_results[[length(merged_all_results) + 1]] <- merged_df
}

# Unir todos los data frames en uno solo
final_combined_table <- bind_rows(merged_all_results)

final_combined_table <- final_combined_table %>%
  left_join(data %>% select(ID, Phylogenetic_Clade), by = c("Strain" = "ID"))

ordered_vars=c("SpB-YPS138-SRR4074411", "SpB-Tank-YMX005537","SpB-Bra-SRR4074412","SpC-SRR7500262",
               "Sp_American","SpB_Undefined","SpA-SRR4074385")

final_combined_table_2 = final_combined_table[,c(1:26,513) ]
colnames(final_combined_table_2)

filtered_data <- final_combined_table %>%
  filter(Phylogenetic_Clade %in% clades_to_plot,
         TrueOrigin %in% ordered_vars)




filtered_data <- final_combined_table %>%
  filter(Phylogenetic_Clade %in% clades_to_plot,
         TrueOrigin %in% ordered_vars,
         !TrueOrigin %in% excluded_variables)



#### jitterplots ####

final_combined_table_2 = final_combined_table[,c(1:26,513) ]

filtered_data <- final_combined_table_2 %>%
  filter(Phylogenetic_Clade %in% clades_to_plot) %>%
  arrange(Phylogenetic_Clade, Strain) %>%
  mutate(
    Strain = factor(Strain, levels = unique(Strain)),  # keep strain order
    Variable = factor(TrueOrigin, levels = ordered_vars) # enforce facet order
  )

write.csv(final_combined_table_2, "Final_combined_table_2.250512.csv")
# All_noOrth de /mnt/Timina/lmorales/Public/ymez/stats/nonMainRefORFs/2024/
# Esto lo busqué en github y me llevó a ArtSace/Fig3/Journal
# Ahí encontré a HetTable_2_240903_conAlp_mod.csv
# De ahí entonces busque en drobpox que es en


final_combined_table_2 = read.csv("Final_combined_table_2.250512.csv")
final_combined_table_2_perGene <- final_combined_table_2 %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = trimws(Genes))  # Remove any leading/trailing whitespace
View(final_combined_table_2_perGene) #Cambiar para que quede por gen.

IntroHet = read.csv("C:/Users/HP/Dropbox/Posdoc/lavis/int/AllSegments/AllIntrogresions/HetTable_2_240903_conAlp_mod.csv")
View(final_combined_table_2)



final_combined_table_2_perGene_Het <- final_combined_table_2_perGene %>%
  left_join(
    IntroHet %>% 
      select(Strain, SAPA_Gene, SACE_Gene, SAPA_cov, SACE_cov, Ratio1, Ratio2),
    by = c("Strain" = "Strain", "Genes" = "SAPA_Gene")
  )
View(final_combined_table_2_perGene_Het) # Agregar los valores de IntroHet

final_combined_table_2_perGene_Het$EsHeterocigoto = final_combined_table_2_perGene_Het$Ratio1 >= 0.25 & final_combined_table_2_perGene_Het$Ratio1 <= 4

final_combined_table_2_perGene_Het
final_combined_table_2
final_combined_table_2_perGene_Het_filtered <- final_combined_table_2_perGene_Het %>%
  filter(TrueIntrogressions == TRUE, TrueOrigin != "Sp_und_outgroup")
final_combined_table_2_filtered <- final_combined_table_2 %>%
  filter(TrueIntrogressions == TRUE, TrueOrigin != "Sp_und_outgroup")

dim(final_combined_table_2)
dim(final_combined_table_2_filtered)
dim(final_combined_table_2_perGene_Het)
dim(final_combined_table_2_perGene_Het_filtered)

#write.csv(final_combined_table_2_perGene, "final_combined_table_2_perGene_250604.csv")
#write.csv(final_combined_table_2_perGene_Het, "final_combined_table_2_perGene_Het_250604.csv")
#write.csv(final_combined_table_2_perGene_Het_filtered, "final_combined_table_2_perGene_Het_filtered_250604.csv")
#write.csv(final_combined_table_2, "Final_combined_table_2_250604.csv")
#write.csv(final_combined_table_2_filtered, "final_combined_table_2_filtered_250604.csv")

# nos dimos cuenta que el All_Blocks estaba incompleto, así que lo volví a cargar y a merge
final_combined_table_2=read.csv("c:/Users/HP/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/figs_code/figsR/Final_combined_table_2_250604.csv")
all_blocks <- read_csv("C:/Users/HP/Desktop/FromIvan/All_blocks_250609.csv")

final_combined_table_2=read.csv("c:/Users/jabra/Dropbox/YMXdomesticated_SACE_PopGen/Draft_Estructura_SACE/figs_code/figsR/Final_combined_table_2_250604.csv")
all_blocks <- read.csv("C:/Users/jabra/Dropbox/Escritorio10072025/FromIvan/All_blocks_250609.csv")


dim(all_blocks)
dim(final_combined_table_2)
merged_genes_expanded <- merged_final_combined_table_2 %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = trimws(Genes))  # Elimina espacios extra si los hay
