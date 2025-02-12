#' Plot a boxplot of a RiboClass object’s counts. 
#' @description This plot is useful to check if the samples are alike in their
#' raw counts.
#' @param ribo A RiboClass object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param outlier Show boxplot outlier values.
#' @param horizontal Show boxplot horizontally.
#' @return A ggplot object.
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' boxplot_count(ribo_toy,"run")
#' 
boxplot_count <- function(ribo, color_col = NA,
                          outlier = TRUE, horizontal = FALSE) {
  ribo_matrix <- extract_data(ribo, "count", position_to_rownames = TRUE)
  
  return(.plot_boxplot_samples(ribo_matrix, "count",
                               ribo[["metadata"]], color_col, outlier,
                               horizontal = horizontal))
}

#' Plot boxplot representing the C-score values of all samples for each
#' individual annotated site.
#' Sites are sorted by their median.
#' @inheritParams boxplot_count 
#' @param sort_by Sort sites by median ("median", default) by variance ("var")
#'  or IQR ("iqr").
#' @return a ggplot geom_boxplot
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' boxplot_cscores(ribo_toy)
#' 
boxplot_cscores <- function(ribo,outlier = TRUE, sort_by = c("median","iqr","var")[1], horizontal = FALSE) { #fonction principale qui appelle la fonction interne .plot_boxplot_sites 
  ribo_m <- extract_data(ribo,only_annotated = TRUE) #extraction des données annotées
  
  if(nrow(ribo_m) == 0) { #Vérification que la riboclass est annotée
    stop("No annotated site found. Please use annotate_site() on your RiboClass before calling this function.")
    }
  
  return(.plot_boxplot_sites(ribo_m, #fonction interne
                      values_to_plot = "cscore",outlier = TRUE, 
                      sort_by = sort_by, # cette variable permet de trier les sites annotés selon la médiane (défaut) ou l'IQR ou la var
                      horizontal = horizontal))
}

#' Internal function for boxplot_cscore
#'
#' @param matrix Sites x Samples C-score/count matrix (output of extract_data()).
#' @param values_to_plot Value to display in plot.
#' @param outlier Show boxplot outlier values.
#' @param horizontal Show boxplot horizontally.
#'
#' @return A ggplot object.
#' @keywords internal
#' Ces deux fonctions au dessus et en dessous permettent de tracer des boxplots pour visualiser la 
#' distribution des C-scores associés aux sites annotés dans un objet RiboClass.


.plot_boxplot_sites <- function(matrix, values_to_plot, outlier, sort_by, horizontal) { #fonction principale appelée dans le Rmd
  site <- cscore <- NULL #init des variables
  id_vars <- "site" #utile pour transformation en format long
    
    #tri des sites selon 3 critères (iqr, median ou var) par défaut: median
    if (tolower(sort_by) == "iqr") { # convertit la variable sort_by en minuscule et vérifie si sort_by = "iqr"
      matrix_melted <- get_IQR(matrix) # si TRUE, la fonction get_IQR est appelée
      #print(head(matrix_melted)) #pour débuguer
      
    } else if(tolower(sort_by) == "median") { #sinon convertit la variable sort_by en minuscule et vérifie si sort_by = "median"
      matrix_melted <- reshape2::melt(matrix, id.vars = id_vars, # Transformation du tableau matrix en format long. id.vars = site
                                      value.name = values_to_plot) # Colonne qui stockera les valeurs ("cscore")
      print(head(matrix_melted))
      
    } else if (tolower(sort_by) == "var")  { # sinon si sort_by = var on tri par var 
      matrix_melted <- get_IQR(matrix, "var") # appel fonction get_IQR 
      #print(head(matrix_melted))
      
    } else {
      stop("Choose either \"median\", \"iqr\" or \"var\" for sort_by param.\n  \"", # message d'erreur si aucun des 3 n'est reconnu 
           sort_by,"\" is not a valid argument.")
    }
  
  shape_outlier <- NA 
  if (outlier) # si outlier = TRUE
    shape_outlier <- 19 # représenté par cercle plein (19)
  
  if ((tolower(sort_by) %in% c("iqr","var"))) { # si sort_by = iqr ou var
    method <- ifelse(tolower(sort_by) =="iqr","IQR","variance") # pour le titre dynamique
   
  p <- ggplot(matrix_melted, aes(x = site, y = cscore)) + # création du plot
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1),
          legend.position = "top",
          axis.title=element_text(size=14),
          axis.text.y = element_text(size=12)) +
    labs(x = paste0("Site (sorted by ",method,")"),
         y = "C-score")
  # Ajout de la courbe des iqr à droite du graphe----------------------------
  # Trier les sites en fonction de l'IQR moyen
  # site_order <- matrix_melted %>%
  #   dplyr::group_by(site) %>%
  #   dplyr::summarize(mean_iqr = mean(iqr, na.rm = TRUE)) %>%
  #   dplyr::arrange(mean_iqr) %>%
  #   dplyr::pull(site)  # Extraire la liste des sites triés
  # 
  # # Appliquer l'ordre des facteurs aux sites
  # matrix_melted$site <- factor(matrix_melted$site, levels = site_order)
  
  # Création du graphique
  # p2 <- ggplot(matrix_melted, aes(x = site, y = iqr, group = 1)) +
  #   geom_point(color = "blue") +     # Points des IQR
  #   geom_line(color = "blue") +      # Ligne reliant les IQR
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  #         axis.title = element_text(size = 14),
  #         axis.text.y = element_text(size = 12)) +
  #   labs(x = "Sites (triés par IQR moyen)", 
  #        y = "IQR", 
  #        title = "Évolution de l'IQR par site")
  p2 <- ggplot(matrix_melted, aes(x = iqr)) +
    geom_density(fill = "blue", alpha = 0.4) +
    coord_flip() +
    theme_minimal() +
    labs(x = "IQR", y = "Density")
  # ------------------------------------------------------------------------------- 
  library(patchwork)
  layout <- c(
    area(t = 0, l = 0, b = 5, r = 4),
    area(t = 2, l = 4, b = 5, r = 5)
  )
   # Empiler les deux graphiques
   final_plot <- p | p2 + # Le `/` permet d'empiler les graphiques verticalement avec {patchwork} et | horizontalement
     plot_layout(design = layout)#plot_layout(widths = c(0.85,0.15)) 
  return(final_plot)
  
  } else { # sinon (si sort_by = median)
    matrix_melted <- reshape2::melt(matrix, id.vars = id_vars, # Trie site par médiane avec stats::reorder()
                                    value.name = values_to_plot)
    
  p <- ggplot2::ggplot(matrix_melted,
                       ggplot2::aes(x = stats::reorder(site, # Trie les valeurs de site en fonction d'une métrique (metric)
                                                       !!rlang::sym(values_to_plot), # récupère dynamiquement la colonne correspondante à values_to_plot ("cscore", "median", "iqr", etc.)
                                                       na.rm = TRUE),
                                    y = !!rlang::sym(values_to_plot))) +
    ggplot2::geom_boxplot(outlier.shape = shape_outlier) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1)) +
    ggplot2::xlab("Site (sorted by median)") +
    ggplot2::ylab("C-score")
  }
  
  if (horizontal) # si horisontal = TRUE
    p <- p + ggplot2::coord_flip() # inverse les axes
  
  return(p)
}

#' (internal) plot boxplot for a given matrix of values
#'
#' @param matrix Sites x Samples C-score/count matrix (output of extract_data()).
#' @param metadata Metadata of samples in matrix.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param outlier Show boxplot outlier values.
#' @param values_col_name Name of the column containing the value (either count
#' or cscore).
#' @param horizontal Show boxplot horizontally.
#'
#' @return ggplot boxplot
#' @keywords internal

.plot_boxplot_samples <- function(matrix, values_col_name, metadata,
                                  color_col = NA, outlier= TRUE, horizontal=TRUE) {
  Sample <- NULL
  id_vars <- "Sample"
  matrix <- log10(matrix) #conversion des c-score en log10
  matrix_inv <- as.data.frame(t(matrix)) #transposer la matrice
  # matrix_inv <- tibble::rownames_to_column(matrix_inv, "Sample")
  matrix_inv["Sample"] <- rownames(matrix_inv)
  
  #vérifier si color_col est fourni ?
  if (!is.na(color_col)) {
    matrix_inv <- cbind(matrix_inv, metadata[color_col])
    id_vars <- c(color_col, "Sample")
  }else{
    # Si color_col est manquant, créer une colonne "qc" pour identifier les outliers
    matrix_inv <- cbind(matrix_inv, qc=0)
    matrix_inv$qc[apply(matrix,2,median,"na.rm"=T)<=2] <- 1 
    id_vars <- c("qc", "Sample")
    
  }
  # Transformer en format "melted" pour ggplot
  matrix_melted <- reshape2::melt(matrix_inv, id.vars = id_vars, #création d'une matrice avec colonnes "qc" (=0 ou 1 si outlier) "Sample", "variable" et "count"
                                  value.name = values_col_name)
  # Transformer "Sample" en facteur pour contrôler l'ordre
  matrix_melted[["Sample"]] <- factor(matrix_melted[["Sample"]],
                                      levels = unique(matrix_melted[["Sample"]]))
  shape_outlier <- NA

  if (outlier)
    shape_outlier <- 19
  
  # Initialiser le plot avec ou sans la couleur définie
  if (is.na(color_col)) {   #s'il n'y a pas de colonne color_col dans les metadata

    p <- ggplot2::ggplot(matrix_melted,
         ggplot2::aes(x = Sample, y = !!sym(values_col_name), #on créer la colonne qc pour mapper les couleurs
                                      fill = factor(qc))) +
         ggplot2::theme(legend.position="none",
         ggplot2::scale_fill_manual(values=c("white", "red")))
  } else {
    p <- ggplot2::ggplot(matrix_melted,
                         ggplot2::aes(x = Sample , y = !!sym(values_col_name), #les couleurs sont mappés à la colonne correspondante
                                      fill = !!sym(color_col)))
  }
  #ajout boxplot
  p <- p + ggplot2::geom_boxplot(outlier.shape = shape_outlier) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1))

  #ajout ligne horizontale
  p <- p + ggplot2::geom_hline(yintercept = 2, colour = "blue")

  #orientation horizontale (optionnel)
  if (horizontal)
    p <- p + ggplot2::coord_flip()

  return(p)
}

# --------------------------------------------------------------------------------  
#test qui fonctionne à supprimer si n-2 partie au dessus fonctionne 
# install.packages("tidyverse")
# install.packages("ggtext")
# library(ggtext)
# library(tidyverse) 
# # Créer les étiquettes avec HTML
# matrix_melted <- matrix_melted %>%
#   mutate(
#     Sample.label = paste0(
#       "<span style='color: ", ifelse(qc == 1, "red", "black"), ";'>",
#       Sample, "</span>"
#     )
#   )
# 
# # Construire le graphique
# p <- ggplot(matrix_melted, aes(x = Sample.label, y = !!sym(values_col_name), fill = factor(qc))) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("white", "red")) +
#   theme(
#     legend.position = "none",
#     axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1) # Permet le rendu HTML
#   )
# 
# # Afficher le graphique
# print(p)
# --------------------------------------------------------------------------------  