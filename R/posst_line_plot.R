#' POSST score line plot
#'
#' @param st Spatial transcriptomics seurat object (after Run create_posst_matrix function).
#' @param niche Spatial niches or clusters.
#' @param score The column name of score, which can be the Neighborhood ("neigh_score"), Aggregation ("agg_score"), and Co-localization Score("coloc_score_...").
#' @param split.by a factor in object metadata to split the plot by (NULL for default).
#' @param split.scales Should scales be fixed("fixed",the default,free("free"),or free in one dimension("free_x","free_y"))?
#' @param nrow number of row to demonstrate when splitting the plot (1 for default).
#' @param ncol number of col to demonstrate when splitting the plot.
#' @param col Vector of colors, each color corresponds to an identity class. A list of 74 colors in RColorBrewer package is pre-saved.
#' @param dot.size.range Range of dot size representing the spots number of each niche. By default the dot.size.range is c(2,8).
#' @param ylim Y axis limits of the plot.
#' @param label.size Text size of labels.
#' @param title Title of the plot.
#' @param title.size Text size of the title.
#'
#' @return The ranking line plot of the POSST score calculated with the POSST package, including the Neighborhood ("neigh_score"), Aggregation ("agg_score"), and Co-localization Score("coloc_score_...").
#' @export
#'
#' @examples plot <- posst_line_plot(st = stobj,
#'                                   niche = 'niche',
#'                                   score = 'neigh_score',
#'                                   split.by = NULL,
#'                                   split.scales = 'fixed',
#'                                   nrow = 1,
#'                                   ncol = 2,
#'                                   col = c('niche1' = 'f60000',...),
#'                                   dot.size.range = c(2,8),
#'                                   ylim = c(0,6),
#'                                   label.size = 2.5,
#'                                   title = 'Neighborhood score ranking',
#'                                   title.size = 12)
#'
posst_line_plot <- function(st = stobj,
                            niche = 'niche',
                            score = 'neigh_score',
                            split.by = NULL,
                            split.scales = 'fixed',
                            nrow = 1,
                            ncol = 2,
                            col = unlist(mapply(brewer.pal,brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors,rownames(brewer.pal.info[brewer.pal.info$category == 'qual',]))),
                            dot.size.range = c(2,8),
                            ylim = c(0,6),
                            label.size = 3,
                            ylabs = 'Mean POSST score',
                            title = 'POSST score',
                            title.size = 12){
  library(Seurat)
  library(tidyverse)
  library(RColorBrewer)
  t_all <- st@meta.data
  t_all$niche <- t_all[[niche]]
  t_all$score <- t_all[[score]]
  t_all <- t_all %>%
  filter(total_neighbor_count == 6) %>%
  group_by(niche) %>%
  mutate(sum=1)%>%
  summarise(sum=sum(sum),
            mean_score=mean(score)) %>%
  ungroup()

t_all <- t_all[order(t_all$mean_score,decreasing = T),]
t_all$niche <- factor(t_all$niche,levels = unique(t_all$niche))

midp <- round(nrow(t_all))/2

t_all[1:midp,'toplab'] <- t_all[1:midp,'niche']
t_all[round(midp+1):nrow(t_all),'botlab'] <- t_all[round(midp+1):nrow(t_all),'niche']

plot <- ggplot(data = t_all,mapping = aes(x = niche, y = mean_score,size = sum,color = niche))+
  geom_point()+
  geom_line(aes(group= 1),size = 1,color = 'black')+
  geom_text(aes(label = toplab), size = label.size, hjust = 0, angle = 45,
            position = position_nudge(x = 0.1, y = 0.02),color = 'black')+
  geom_text(aes(label = botlab), size = label.size, hjust = 1, angle = 45,
            position = position_nudge(x = -0.1, y = -0.02),color = 'black')+
  labs(x = '',y = ylabs,size = 'Spots number',title = title)+
  ylim(ylim)+
  scale_color_manual(values = col)+
  scale_size_continuous(range = dot.size.range)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size = title.size))+
  guides(color = 'none')
if (is.null(split.by)){
  plot <- plot
}else{
  t_all <- st@meta.data
  t_all$niche <- t_all[[niche]]
  t_all$score <- t_all[[score]]
  t_all$split.by <- t_all[[split.by]]
  t_all <- t_all %>%
    filter(total_neighbor_count == 6) %>%
    group_by(niche,split.by) %>%
    mutate(sum=1)%>%
    summarise(sum=sum(sum),
              mean_score=mean(score)) %>%
    ungroup()

  t_all <- t_all[order(t_all$mean_score,decreasing = T),]
  t_all$niche <- factor(t_all$niche,levels = unique(t_all$niche))
  midp <- round(nrow(t_all))/2

  t_all[1:midp,'toplab'] <- t_all[1:midp,'niche']
  t_all[round(midp+1):nrow(t_all),'botlab'] <- t_all[round(midp+1):nrow(t_all),'niche']

  plot <- ggplot(data = t_all,mapping = aes(x = niche, y = mean_score,size = sum,color = niche))+
    geom_point()+
    geom_line(aes(group= 1),size = 1,color = 'black')+
    geom_text(aes(label = toplab), size = label.size, hjust = 0, angle = 45,
              position = position_nudge(x = 0.1, y = 0.02),color = 'black')+
    geom_text(aes(label = botlab), size = label.size, hjust = 1, angle = 45,
              position = position_nudge(x = -0.1, y = -0.02),color = 'black')+
    labs(x = '',y = ylabs,size = 'Spots number',title = title)+
    ylim(ylim)+
    #scale_color_manual(values = col)+
    scale_size_continuous(range = dot.size.range)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
    theme(plot.title = element_text(hjust = 0.5,size = title.size))+
    guides(color = 'none')+
    facet_wrap(~split.by,scales = split.scales, nrow = nrow, ncol = ncol)
}
return(plot)
}

