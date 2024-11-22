% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posst_line_plot.R
\name{posst_line_plot}
\alias{posst_line_plot}
\title{POSST score line plot}
\usage{
posst_line_plot(
  st = stobj,
  niche = "niche",
  score = "neigh_score",
  split.by = NULL,
  split.scales = "fixed",
  nrow = 1,
  ncol = 2,
  col = unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == "qual",
    ]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == "qual", ]))),
  dot.size.range = c(2, 8),
  ylim = c(0, 6),
  label.size = 3,
  ylabs = "Mean POSST score",
  title = "POSST score",
  title.size = 12
)
}
\arguments{
\item{st}{Spatial transcriptomics seurat object (after Run create_posst_matrix function).}

\item{niche}{Spatial niches or clusters.}

\item{score}{The column name of score, which can be the Neighborhood ("neigh_score"), Aggregation ("agg_score"), and Co-localization Score("coloc_score_...").}

\item{split.by}{a factor in object metadata to split the plot by (NULL for default).}

\item{split.scales}{Should scales be fixed("fixed",the default,free("free"),or free in one dimension("free_x","free_y"))?}

\item{nrow}{number of row to demonstrate when splitting the plot (1 for default).}

\item{ncol}{number of col to demonstrate when splitting the plot.}

\item{col}{Vector of colors, each color corresponds to an identity class. A list of 74 colors in RColorBrewer package is pre-saved.}

\item{dot.size.range}{Range of dot size representing the spots number of each niche. By default the dot.size.range is c(2,8).}

\item{ylim}{Y axis limits of the plot.}

\item{label.size}{Text size of labels.}

\item{title}{Title of the plot.}

\item{title.size}{Text size of the title.}
}
\value{
The ranking line plot of the POSST score calculated with the POSST package, including the Neighborhood ("neigh_score"), Aggregation ("agg_score"), and Co-localization Score("coloc_score_...").
}
\description{
POSST score line plot
}
\examples{
plot <- posst_line_plot(st = stobj,
                                  niche = 'niche',
                                  score = 'neigh_score',
                                  split.by = NULL,
                                  split.scales = 'fixed',
                                  nrow = 1,
                                  ncol = 2,
                                  col = c('niche1' = 'f60000',...),
                                  dot.size.range = c(2,8),
                                  ylim = c(0,6),
                                  label.size = 2.5,
                                  title = 'Neighborhood score ranking',
                                  title.size = 12)

}
