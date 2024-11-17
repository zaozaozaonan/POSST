#' Co-localization score heatmap of all the niches
#'
#' @param st Spatial transcriptomics seurat object (after Run create_posst_matrix function).
#' @param niche Spatial niches or clusters.
#' @param subgroup a factor in object metadata to subset the stobj (e.g. the treatment options, or the sampling sites. NULL for default).
#' @param subgroup_type the a value of subgroup to choose (e.g. treat/control, or tumor/para-tumor. Single value or a vector of values can be set as the input. NULL for default).
#' @param col colors used to plot (rev(brewer.pal(n = 7, name ="RdYlGn") for default).
#' @param scale character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row","column", and "none".
#' @param title Title of the plot.
#'
#' @return The Co-localization heatmap of all the niches, calculating with the POSST package.
#' @export
#'
#' @examples plot <- posst_coloc_heatmap(st = stobj,
#'                                       niche = 'niche',
#'                                       subgroup = 'sampling_site',
#'                                       subgroup_type = c('T','PT'),
#'                                       col = rev(brewer.pal(n = 7, name ="RdYlGn")),
#'                                       scale = 'none',
#'                                       title = 'Niche co-localization score')
posst_coloc_heatmap <- function(st = stobj,
                                niche = 'niche',
                                subgroup = NULL,
                                subgroup_type = NULL,
                                col = rev(brewer.pal(n = 7, name ="RdYlGn")),
                                scale = 'none',
                                title = 'niche co-localization score'
                                ){
  library(Seurat)
  library(tidyverse)
  library(RColorBrewer)
  t_all <- st@meta.data
  t_all$niche <- t_all[[niche]]
  if (is.null(subgroup)){
    t_all <- t_all %>%
      filter(total_neighbor_count == 6) %>%
      group_by(niche) %>%
      select(all_of(paste0('coloc_score_',unique(t_all$niche)))) %>%
      summarise_all(mean) %>%
      ungroup() %>%
      column_to_rownames(.,var = 'niche')
  }else{
    t_all$subgroup <- t_all[[subgroup]]
    t_all <- t_all %>%
      filter(total_neighbor_count == 6) %>%
      filter(subgroup %in% subgroup_type) %>%
      group_by(niche) %>%
      select(all_of(paste0('coloc_score_',unique(t_all$niche)))) %>%
      summarise_all(mean) %>%
      ungroup() %>%
      column_to_rownames(.,var = 'niche')
  }
  plot <- pheatmap::pheatmap(t_all,main = title,cellwidth = 12,cellheight = 12,angle_col = 315,scale = scale,
                     color = colorRampPalette(col)(100))
  return(plot)
}
