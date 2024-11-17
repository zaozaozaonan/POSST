#' Co-localization score lollipop plot between the center niche and other niches
#'
#' @param st Spatial transcriptomics seurat object (after Run create_posst_matrix function).
#' @param niche Spatial niches or clusters.
#' @param center_niche The center niche.
#' @param subgroup a factor in object metadata to subset the stobj (e.g. the treatment options, or the sampling sites. NULL for default).
#' @param subgroup_type the a value of subgroup to choose (e.g. treat/control, or tumor/para-tumor. Single value or a vector of values can be set as the input. NULL for default).
#' @param col colors used to plot (NULL for default).
#' @param dot.size dot size of the lollipop. By default the dot.size is 5.
#' @param ylim Y axis limits of the plot.
#' @param label.size Text size of labels.
#' @param title Title of the plot.
#' @param title.size Text size of the title.
#'
#' @return The lollipop plot of co-localization score between the center niche X and other niches calculated with the POSST package.
#' @export
#'
#' @examples plot <- posst_coloc_lollipop_plot(st = stobj,
#'                                             niche = 'niche',
#'                                             center_niche = 'niche_1',
#'                                             subgroup = 'sampling_site',
#'                                             subgroup = c('T','PT'),
#'                                             col = c('niche1' = 'f60000',...),
#'                                             dot.size = 5,
#'                                             ylim = c(0,6),
#'                                             label.size = 3,
#'                                             title = paste('Co-localization score between',center_niche,'and other niches',subgroup_type,sep = ' '),
#'                                             title.size = 12)
posst_coloc_lollipop_plot <- function(st = stobj,
                                      niche = 'niche',
                                      center_niche = 'niche_1',
                                      subgroup = NULL,
                                      subgroup_type = NULL,
                                      col = NULL,
                                      dot.size = 5,
                                      ylim = c(0,6),
                                      label.size = 3,
                                      title = paste('Co-localization score between',center_niche,'and other niches',subgroup_type,sep = ' '),
                                      title.size = 12){
  library(Seurat)
  library(tidyverse)
  library(RColorBrewer)
  t_all <- st@meta.data
  t_all$niche <- t_all[[niche]]
  if (is.null(subgroup)){
    t_all <- t_all %>%
      filter(total_neighbor_count == 6) %>%
      filter(niche == center_niche) %>%
      select(all_of(paste0('coloc_score_',unique(t_all$niche)))) %>%
      summarise_all(mean)
      }else{
    t_all$subgroup <- t_all[[subgroup]]
    t_all <- t_all %>%
      filter(total_neighbor_count == 6) %>%
      filter(niche == center_niche) %>%
      filter(subgroup %in% subgroup_type) %>%
      select(all_of(paste0('coloc_score_',unique(t_all$niche)))) %>%
      summarise_all(mean)
    }
  t_all <- data.frame(t(t_all))
  colnames(t_all) <- 'mean_coloc_score_'
  t_all <- rownames_to_column(t_all,var = 'coloc_niche')
  t_all <- t_all[order(t_all$mean_coloc_score_,decreasing = T),]
  t_all$coloc_niche <- factor(t_all$coloc_niche,levels = t_all$coloc_niche,labels = gsub('coloc_score_','',t_all$coloc_niche))

  plot <- ggplot(data = t_all,mapping = aes(x =  coloc_niche, y = mean_coloc_score_ ,color = coloc_niche))+
    geom_point(size = dot.size)+
    geom_segment(aes(x = coloc_niche ,xend = coloc_niche ,y =  0,yend = mean_coloc_score_ ))+
    geom_text(aes(label = coloc_niche), size = label.size, hjust = 0, angle = 90,
              position = position_nudge(x = 0, y = 0.15),color = 'black')+
    labs(x = '',y = 'Mean co-localization score',title = title)+
    ylim(ylim)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size = title.size))+
    theme(panel.grid = element_blank())+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+ guides(color = 'none')

  if (is.null(col)){
    plot <- plot
  }else{
    plot <- plot +
      scale_color_manual(values = col)
  }
  return(plot)
}

