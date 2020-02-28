#' DCLEAR
#'
#' @import Matrix
#' @importFrom matrixStats rowSds rowVars rowMedians rowMins rowMaxs
#' @import futile.logger
#' @importFrom phangorn NJ treedist RF.dist
#' @importFrom igraph as.igraph distances graph.adjacency ends incident_edges adjacent_vertices degree permute.vertices vcount V contract simplify set_vertex_attr graph.empty add_edges add_vertices
#' @importFrom purrr map
#' @importFrom stringr str_pad
#' @importFrom abind abind
#' @importFrom Biostrings BStringSet oligonucleotideFrequency
#' @importFrom dplyr filter sample_n mutate select collect
#'
NULL

#' register S3 classes
setOldClass('phyDat')
setOldClass('phylo')
setOldClass('igraph')

#' register new S4 classes
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
		dropout = 'logical'
	)
)
