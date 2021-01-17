#' DCLEAR: A package for DCLEAR: Distance based Cell LinEAge Reconstruction
#'
#' Distance based methods for inferring lineage treess from single cell data
#' 
#' @import tidyverse
#' @import methods
#' @import Matrix
#' @import dplyr
#' @importFrom matrixStats rowLogSumExps
#' @importFrom igraph distances graph.adjacency degree permute.vertices vcount V contract simplify set_vertex_attr
#' @importFrom purrr map
#' @importFrom stringr str_pad
#' @importFrom BiocParallel bplapply
#' @importFrom phangorn phyDat
#' @importFrom stats dist
#' @importFrom Rcpp evalCpp
#' @importFrom ape write.tree read.tree
#' @importFrom tidyr replace_na
#' @importFrom keras keras_model_custom
#' @useDynLib DCLEAR
#' @docType package
#' @name DCLEAR
#'
NULL
# > NULL

