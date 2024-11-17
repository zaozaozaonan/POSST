#' Find the boundary between the center niche and other niches
#'
#' @param st Spatial transcriptomics seurat object.
#' @param niche Spatial niches or clusters.
#' @param center_niche The center niche.
#' @param images Vector of the name for all the images. By default, all the slides are included.
#' @param strictness Number of neighboring spots in boundary finding, indicating the strictness of the boundary definition. This parameter can take the values 1, 2, or 3, with 3 representing the strictest setting. (By default strictness = 2)
#' @param platform Platform of spatial transcriptomics. 10x visium or BMK s1000 (L13 bin) spatial transcriptomics can be loaded as input for this package. Corresponding values are "10x" and "s1000". By default, platform = '10x'.
#'
#' @return The boundary spots between the center niche and other niches were labeled in "boundary". The minimal spatial distance for each spots to the boundary ("MDB" for short) were calculated. "MDBu" represents the value adjusted from MDB, which holds practical significance and is measured in units of 100 micrometers (μm).
#' @export
#'
#' @examples st <- find_posst_boundary(st = stobj,
#'                                     niche = 'niche',
#'                                     center_niche = 'niche_1',
#'                                     images = names(st@images),
#'                                     strictness = 2,
#'                                     platform = '10x')
find_posst_boundary <- function(st = stobj,
                          niche = 'niche',
                          center_niche = 'niche_1',
                          images = names(st@images),
                          strictness = 2,
                          platform = '10x'){
  library(Seurat)
  library(tidyverse)
  m_niche <- st[[niche]];colnames(m_niche) = 'niche'
  m_niche$ts <- ifelse(m_niche$niche == center_niche,'center','other')

  for (imm in 1:length(images)){
    co <- st@images[[images[imm]]]@coordinates
    m_niche_s <- m_niche[rownames(co),]
    co <- merge(co,m_niche_s,by = 'row.names')
    if(platform=='10x'){
      for(i in 1:nrow(co)){
        x <- co[i,3]
        y <- co[i,4]
        n <- co%>%
          filter(row%in%c(x-1,x,x+1))%>%
          filter(col%in%c(y-2,y-1,y+1,y+2))#filter the surrending matrix for each spot
        coc <- nrow(n%>%filter(ts=="center"))
        coo <- nrow(n%>%filter(ts=="other"))
        co[i,'boundary'] <- ifelse(coc>=strictness&coo>=strictness,'BD',co[i,'niche'])
      }
    }else if(platform=='s1000'){
      co$row <- round(co$row/0.6073109)
      co$col <- round(co$col/0.3506311)
      for(i in 1:nrow(co)){
        x <- co[i,3]
        y <- co[i,4]
        n <- co%>%
          filter(row%in%c(x-25,x-13,x-12,x+12,x+13,x+25))%>%
          filter(col%in%c(y-38,y-37,y-1,y+1,y+37,y+38))#filter the surrending matrix for each spot
        coc <- nrow(n%>%filter(ts=="center"))
        coo <- nrow(n%>%filter(ts=="other"))
        co[i,'boundary'] <- ifelse(coc>=strictness&coo>=strictness,'BD',co[i,'niche'])
      }
    }
    #计算MDB and MDBu
    for(a in 1:nrow(co)){
      x <- co[a,5]
      y <- co[a,6]
      ed <- function(xy){
        x1 <- xy[1]
        y1 <- xy[2]
        sqrt((x-x1)^2+(y-y1)^2)
      }
      bdm <- co%>%filter(grepl("BD",x = co$boundary))
      dtb <- apply(bdm[,c(5,6)],1,ed)
      mdb <- min(dtb)
      co[a,"MDB"] <- ifelse(co[a,"ts"]=='center',-mdb,mdb)
    }
    co['maxmdb_abs'] <- max(abs(co$MDB))
    co['minmdb_abs'] <- min(abs(co %>% filter(MDB!=0) %>% select(MDB)))
    co['MDBu'] <- round(co$MDB/abs(co$minmdb_abs),digits = 2)
    if(imm==1){b <- co}else{b <- rbind(b,co)}
  }
  t_all <- st@meta.data
  t_all['boundary'] <- plyr::mapvalues(rownames(t_all),from = b$Row.names,to = b$boundary)
  t_all['MDB']<- plyr::mapvalues(rownames(t_all),from = b$Row.names,to = b$MDB) %>% as.numeric()
  t_all['MDBu']<- plyr::mapvalues(rownames(t_all),from = b$Row.names,to = b$MDBu) %>% as.numeric()
  st@meta.data <- t_all
  return(st)
}
