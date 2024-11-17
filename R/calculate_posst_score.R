#' Create the Matrix for Calculating Neighborhood/Aggregation/Co-localization Score
#'
#' @param st Spatial transcriptomics seurat object.
#' @param niche Spatial niches or clusters.
#' @param images Vector of the name for all the images. By default, all the slides are included.
#' @param platform Platform of spatial transcriptomics. 10x visium or BMK s1000 (L13 bin) spatial transcriptomics can be loaded as input for this package. Corresponding values are "10x" and "s1000". By default, platform = '10x'.
#'
#' @return The total number of spots around the center spot and neighborhood/aggregation/co-localization score. Corresponding values are "total_neighbor_count", "neigh_score", "agg_score", "coloc_score_..."
#' @export
#'
#' @examples
#' images <- names(stobj@images)
#' st <- calculate_posst_score(st = stobj,
#'                             niche = 'niche',
#'                             images = images,
#'                             platform = 10x)
calculate_posst_score <- function(st = stobj,
                                  niche = 'niche',
                                  images = names(st@images),
                                  platform = '10x') {
  library(Seurat)
  library(tidyverse)
  m_niche <- st[[niche]];colnames(m_niche) = 'niche'

  for (imm in 1:length(images)){
    co <- st@images[[images[imm]]]@coordinates#get coordincates matrix
    m_niche_s <- m_niche[rownames(co),]
    co <- merge(co,m_niche_s,by = 'row.names')#add niche information to the coor matrix
    if(platform == '10x'){
      for(i in 1:nrow(co)){
        x <- co[i,'row']
        y <- co[i,'col']
        n <- co%>%
          filter(row%in%c(x-1,x,x+1))%>%
          filter(col%in%c(y-2,y-1,y+1,y+2))#filter the surrending matrix for each spot
        co[i,'total_neighbor_count'] <- nrow(n)
        co[i,'neigh_score'] <- length(unique(n$niche))
        co[i,'agg_score'] <- nrow(n %>% filter(niche == co[i,'niche']))
        for (nic in levels(co$niche)){
          co[i,paste0('coloc_score_',nic)] <- nrow(n %>% filter(niche == nic))
        }
      }
    }else if (platform == 's1000'){
      # mindist <- dist(cor$row)
      # mindist <- as.matrix(mindist)
      # mindist <- mindist[mindist!=0]
      # min(mindist)#0.6073109
      co$row <- round(co$row/0.6073109)
      co$col <- round(co$col/0.3506311)
      for(i in 1:nrow(co)){
        x <- co[i,'row']
        y <- co[i,'col']
        n <- co%>%
          filter(row%in%c(x-25,x-13,x-12,x+12,x+13,x+25))%>%
          filter(col%in%c(y-38,y-37,y-1,y+1,y+37,y+38))#filter the surrending matrix for each spot
        co[i,'total_neighbor_count'] <- nrow(n)
        co[i,'neigh_score'] <- length(unique(n$niche))
        co[i,'agg_score'] <- nrow(n %>% filter(niche == co[i,'niche']))
        for (nic in levels(co$niche)){
          co[i,paste0('coloc_score_',nic)] <- nrow(n %>% filter(niche == nic))
        }
      }
    }
    if(imm==1){nbor <- co}else{nbor <- rbind(nbor,co)}
    st@meta.data[,colnames(nbor)[8:ncol(nbor)]] <- nbor[,8:ncol(nbor)]
    return(st)
  }
}

