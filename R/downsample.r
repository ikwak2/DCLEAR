#' downsample
#'
#' Sample a lineage tree
#'
#' @param x a lineage_tree object
#' @param n number of leaves (tips) in the down-sampled tree
#' @param ... additional parameters
#'
#' @return a lineage_tree object
#'
#' @export
#'
setMethod(
	'downsample',
	signature(
		x = 'lineage_tree'
	),
	function(
		x,
		n = 10L,
		...
	){

		is_leaf <- degree(x@graph, mode = 'out') == 0
		if (n < sum(is_leaf)){
	  	leaf <- sample(names(x@x)[is_leaf], 1L)
			d <- distances(x@graph, v = leaf, to = names(which(is_leaf)))
			leaves <- colnames(d)[order(d)[1:n]]
#			leaves <- sample(colnames(d), n)
			x <- subtree(x, leaves)
		}
		x
	}
) # downsample


#' downsample
#'
#' Sample a lineage tree
#'
#' @param x a igraph object
#' @param n number of leaves (tips) in the down-sampled tree
#' @param ... additional parameters
#'
#' @return a phylo object
#'
#' @export
#'
setMethod(
	'downsample',
	signature(
		x = 'igraph'
	),
	function(
		x,
		n = 10L,
		...
	){

		browser()

		is_leaf <- degree(x, mode = 'out') == 0
		if (n <= sum(is_leaf)){
	  	leaves <- sample(unlist(V(x)$name[is_leaf]), n)
			d <- distances(x, leaves, mode = 'in')
		  v <- V(x)$name[(!is.infinite(d)) %>% colSums() > 0]  # subgraphs that connect to the leaves
			v <- unlist(v)
			x <- induced_subgraph(x, v)
		}
		x
	}
)

