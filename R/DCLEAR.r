#' DCLEAR: A package for DCLEAR: Distance based Cell LinEAge Reconstruction
#'
#' Distance based methods for inferring lineage treess from single cell data
#' 
#' @import tidyverse
#' @import dplyr 
#' @import Matrix
#' @importFrom matrixStats rowLogSumExps
#' @import futile.logger
#' @importFrom igraph distances graph.adjacency degree permute.vertices vcount V contract simplify set_vertex_attr
#' @importFrom purrr map
#' @importFrom stringr str_pad
#' @importFrom BiocParallel bplapply
#' @importFrom phangorn phyDat
#' @importFrom stats dist
#' @importFrom Rcpp evalCpp
#' @importFrom ape write.tree read.tree
#' @importFrom tidyr replace_na
#' @useDynLib DCLEAR
#' @docType package
#' @name DCLEAR
#'
NULL
# > NULL

