#' DCLEAR: A package for DCLEAR: Distance based Cell LinEAge Reconstruction
#'
#' Distance based methods for inferring lineage treess from single cell data
#' 
#' @import tidyverse
#' @import methods
#' @import Matrix
#' @import dplyr
#' @importFrom matrixStats rowLogSumExps
#' @importFrom igraph distances graph.adjacency degree permute.vertices vcount V V<- E E<- contract simplify set_vertex_attr graph_from_edgelist induced_subgraph topo_sort
#' @importFrom purrr map
#' @importFrom stringr str_pad
#' @importFrom BiocParallel bplapply
#' @importFrom phangorn phyDat RF.dist
#' @importFrom stats dist as.dist rgamma runif
#' @importFrom Rcpp evalCpp
#' @importFrom ape write.tree read.tree fastme.bal
#' @importFrom tidyr replace_na
#' @importFrom rlang .data
#' @importFrom rBayesianOptimization BayesianOptimization
#' @importFrom BiocGenerics rbind
#' @useDynLib DCLEAR
#' @docType package
#' @name DCLEAR
#'
NULL
# > NULL

