#' DCLEAR: A package for DCLEAR: Distance based Cell LinEAge Reconstruction
#'
#' Distance based methods for inferring lineage treess from single cell data
#' 
#' @import tidyverse
#' @import dplyr 
#' @import Matrix
#' @importFrom matrixStats rowSds rowVars rowMedians rowMins rowMaxs rowProds logSumExp
#' @import futile.logger
#' @importFrom phangorn NJ treedist RF.dist
#' @importFrom igraph as.igraph distances graph.adjacency ends incident_edges adjacent_vertices degree permute.vertices vcount V contract simplify set_vertex_attr graph.empty add_edges add_vertices
#' @importFrom purrr map
#' @importFrom stringr str_pad
#' @importFrom Biostrings BStringSet oligonucleotideFrequency
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

setOldClass('phyDat')
setOldClass('phylo')
setOldClass('igraph')

setClass(
	'lineage_tree',
	representation(
		x = 'phyDat',
		graph = 'igraph',
		outcome_prob = 'numeric',
		alphabets = 'character',
		division = 'numeric',
		n_samples = 'numeric',
		n_targets = 'numeric',
		deletion = 'logical',
		dropout_prob = 'numeric'
	)
)

setClass(
	'kmer_summary',
	representation(
		df = 'data.frame',
    alphabets = 'character',
    kmers = 'character',
    outcome_prob = 'numeric',
    sequence_length = 'numeric',
    division = 'numeric',
    reps = 'numeric',
    max_distance = 'numeric',
    mutation_prob = 'numeric',
    k = 'numeric',
		dropout_prob = 'numeric'
	)
)
