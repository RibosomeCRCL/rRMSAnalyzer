# ==================================================

#' Plot a correlation heatmap from a riboclass object.
#'  
#' Shows the correlation **distance** between samples.
#' @md
#' @inheritParams plot_heatmap
#' @param values_col Name of the column containing the value (either count or cscore).
#' @return ComplexHeatmap object
#' @export
#'
#' @examples 
#' data("ribo_toy")
#' plot_heatmap_corr(ribo_toy,"count","condition")
#' 
plot_heatmap_corr <- function(ribo, values_col, color_col=NULL) {
  matrix <- extract_data(ribo, values_col, position_to_rownames = TRUE)
  if(!is.null(color_col)) check_metadata(ribo,color_col)
  .plot_heatmap_corr(matrix, ribo[["metadata"]], color_col = color_col)
}
#' Internal function of plot_heatmap_corr.
#'
#' @param cscore_matrix  Sites x Samples C-score matrix (output of extract_data()).
#' @param metadata Metadata of samples in matrix
#' @param color_col Vector of the metadata columns’ name used for coloring samples.
#'
#' @return ComplexHeatmap heatmap
#' @keywords internal
#' 
.plot_heatmap_corr <- function(cscore_matrix, metadata,
                               color_col) {
  
  if(!is.null(color_col)) { # si aucune couleur est donnée
    col <- generate_palette(metadata,color_col) # génère une palette de couleur pour le heatmap
    column_ha <- ComplexHeatmap::HeatmapAnnotation(df = metadata[color_col], col = col) #créer les annotations de la heatmap
  } else { # sinon
    column_ha <- NULL
  }
    corr_matrix <- stats::cor(cscore_matrix,use = "complete.obs") # calcule la corrélation entre les colonnes ou lignes de la matrice, la corrélation est calculée en ignorant les valeurs manquantes (NA)
    pearson_color <- colorRamp2::colorRamp2(c(0,0.5,1),c("red","white", "blue")) # 0 = rouge, 0.5 = blanc, 1 = bleu
   
 ht <- ComplexHeatmap::Heatmap(corr_matrix,col = pearson_color,name = "Pearson correlation",
                          row_title = "Sample",column_title = "Pearson correlation between samples", 
                          column_title_side = "top", #Position du titre des colonnes en haut
                          cluster_rows = FALSE, cluster_columns = TRUE, # Pas de clustering des lignes (les échantillons restent dans l'ordre original), Clustering des colonnes (les échantillons corrélés seront regroupés)
                          clustering_distance_columns = "pearson", # Utilisation de la distance de Pearson pour mesurer les similarités entre échantillons
                          heatmap_legend_param = list(
                            legend_direction = "horizontal" ),
                          clustering_method_columns = "ward.D2", # Utilisation de la méthode de clustering Ward.D2, qui regroupe les échantillons minimisant la variance intra-groupe
                          column_split = 3, # Divise les colonnes en 3 groupes distincts après clustering
                          top_annotation = column_ha, # ajoute les annotations définis dans la boucle if plus haut
                          row_names_gp = grid::gpar(fontsize = 6))
  
 ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom") # Affiche la légende horizontalement
  
}
