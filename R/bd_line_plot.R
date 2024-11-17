#' Visualize spatial clustering and expression data distribution associated with the boundary.
#'
#'
#' @param st Spatial transcriptomics seurat object.
#' @param group.by The categorical or numeric features for the spatial distribution visualization.
#' @param features Multiple numeric features for the spatial distribution visualization.
#' @param dis.group The distance group. A vector to split the spatial distance associated with boundary (e.g. c(-Inf,-4,-3,-2,-1,0,1,2,3,4,Inf)).
#' @param sample_id If multiple sample is merged in the st object, the column of different sample_id should be filled in here. By default, sample_id = 'orig.ident'.
#' @param subgroup If only a subset of the st object is used to plot, fill in the column name of the categorical variable here and the spots keep in for "subgroup_type".
#' @param subgroup_type See subgroup (NULL for default).
#' @param ylim Y axis limits of the plot (NULL for default).
#' @param method Use the mean or median of the features value among spots in certain distance group. By default, method = "median". "mean" can also be used when the features conforms to a normal distribution.
#' @param col colors used to plot (NULL for default).
#' @param title Title of the plot.
#' @param title.size Text size of the title.
#' @param slot If plotting a feature, which data slot to pull from (counts, data, or scale.data).
#' @param min.cutoff Vector of minimum and maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param max.cutoff See min.cutoff
#' @param shape Control the shape of the spots - same as the ggplot2 parameter. The default is 21, which plots circles - use 22 to plot squares.
#' @param pt.size Size of the spots. The default is 0.
#'
#' @return Line plot of features or discrete grouping (e.g. cluster assignments) distribution associated with the boundary.
#' @export
#'
#' @examples plot <- bd_line_plot(st = stobj,
#'                                group.by = NULL,
#'                                features = NULL,
#'                                slot = 'data',
#'                                min.cutoff = NA,
#'                                max.cutoff = NA,
#'                                dis.group = c(-Inf,-3,-2,-1,0,1,2,3,Inf),
#'                                sample_id = 'orig.ident',
#'                                subgroup = NULL,
#'                                subgroup_type = NULL,
#'                                ylim = NULL,
#'                                method = 'median',
#'                                col = NULL,
#'                                shape = 21,
#'                                pt.size = 0,
#'                                title = 'Spatial distribution',
#'                                title.size = 12)

bd_line_plot <- function(st = stobj,
                         group.by = NULL,
                         features = NULL,
                         slot = 'data',
                         min.cutoff = NA,
                         max.cutoff = NA,
                         dis.group = c(-Inf,-3,-2,-1,0,1,2,3,Inf),
                         sample_id = 'orig.ident',
                         subgroup = NULL,
                         subgroup_type = NULL,
                         ylim = NULL,
                         method = 'median',
                         col = NULL,
                         shape = 21,
                         pt.size = 0,
                         title = 'Spatial distribution',
                         title.size = 12
){
  library(Seurat)
  library(tidyverse)
  library(RColorBrewer)
  t_all <- st@meta.data
  t_all$sample_id <- t_all[[sample_id]]
  t_all$MDBu_c <- cut(t_all$MDBu,breaks = dis.group,include.lowest = TRUE)


  if (is.null(subgroup)){
    t_all <- t_all

  }else{
    t_all$subgroup <- t_all[[subgroup]]
    t_all <- t_all %>% filter(subgroup%in%subgroup_type)

  }

  cells <- row.names(t_all)

  if(is.null(features)){
    t_all$group.by <- t_all[[group.by]]
      t_all_table <- data.frame(table(t_all$group.by,t_all$sample_id,t_all$MDBu_c))
      colnames(t_all_table) <- c('group.by','sample_id','MDBu_c','Freq')
      t_all_table <- t_all_table %>%
        group_by(sample_id,MDBu_c) %>%
        mutate(percent = Freq/sum(Freq)) %>%
        ungroup()
      t_all_table[grep('NaN',t_all_table$percent),'percent'] <- NA
      if(method == 'mean'){
        t_all_table <- t_all_table %>%
          group_by(group.by,MDBu_c) %>%
          mutate (avg_score = mean(na.omit(percent))) %>%
          ungroup()
      }else if(method == 'median'){
        t_all_table <- t_all_table %>%
          group_by(group.by,MDBu_c) %>%
          mutate (avg_score = median(na.omit(percent))) %>%
          ungroup()
      }

      plot <- ggplot(data = t_all_table,mapping = aes(x = MDBu_c,y = percent,color = group.by))+
        geom_point(data = t_all_table,mapping = aes(x = MDBu_c,y = percent,color = group.by),shape = shape,size = pt.size)+
        geom_point(data = t_all_table,mapping = aes(x = MDBu_c, y = avg_score,color = group.by),shape = 5) +
        geom_line(data = t_all_table,mapping = aes(x = MDBu_c, y = avg_score,group = group.by,color = group.by),size = 1)+
        theme_bw()+
        theme(panel.grid = element_blank())+
        labs(x = '',y = 'Percent',color = group.by)
  }else{
    #Features
    data <- FetchData(object = st, vars = features, cells = cells,
                      layer = slot, clean = FALSE)
    features <- colnames(x = data)

    min.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = min(data[,
                                                             feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = max(data[,
                                                             feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features,
                                                min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }

    data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index],
                             data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index],
                             data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      return(data.feature)
    })
    colnames(x = data) <- features
    rownames(x = data) <- cells
    t_all <- t_all[,!colnames(t_all) %in% features]

    t_all <- cbind(t_all,data)

    if(method == 'mean'){
      t_all <- t_all %>%
        select(sample_id,MDBu_c,all_of(features)) %>%
        reshape2::melt(.,
                       id.vars = c('sample_id','MDBu_c'),
                       measure.vars = features,
                       variable.name = 'features',
                       value.name = 'value') %>%
        group_by(features,sample_id,MDBu_c) %>%
        summarise(avg_score_persam = mean(value)) %>%
        ungroup()%>%
        group_by(features,MDBu_c) %>%
        mutate (avg_score = mean(na.omit(avg_score_persam))) %>%
        ungroup()
    }else if(method == 'median'){
      t_all <- t_all %>%
        select(sample_id,MDBu_c,all_of(features)) %>%
        reshape2::melt(.,
                       id.vars = c('sample_id','MDBu_c'),
                       measure.vars = features,
                       variable.name = 'features',
                       value.name = 'value') %>%
        group_by(features,sample_id,MDBu_c) %>%
        summarise(avg_score_persam = median(value)) %>%
        ungroup()%>%
        group_by(features,MDBu_c) %>%
        mutate (avg_score = median(na.omit(avg_score_persam))) %>%
        ungroup()
    }
    plot <- ggplot(data = t_all,mapping = aes(x = MDBu_c,y = avg_score,color = features))+
      geom_point(data = t_all,mapping = aes(x = MDBu_c,y = avg_score_persam,color = group.by),shape = shape,size = pt.size)+
      geom_point(data = t_all,mapping = aes(x = MDBu_c, y = avg_score,color = features),shape = 5) +
      geom_line(data = t_all,mapping = aes(x = MDBu_c, y = avg_score,group = features,color = features),size = 1)+
      theme_bw()+
      theme(panel.grid = element_blank())+
      labs(x = '',y = method)
  }

  if(is.null(col)){
    plot <- plot
  }else{
    plot <- plot+
      scale_color_manual(values = col)
  }

  if(is.null(ylim)){
    plot <- plot
  }else{
    plot <- plot+
      ylim(ylim)
  }

  if(is.null(title)){
    plot <- plot
  }else{
    plot <- plot+
      labs(title = title)+
      theme(plot.title = element_text(hjust = 0.5,size = title.size))
  }

  return(plot)
}
