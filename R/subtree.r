#' subtree
#'
#' Extract a subtree with specific leaves
#'
#' @param x a lineage_tree object
#' @param leaves leaves of the extracted tree
#' @param ... additional parameters
#'
#' @return a lineage_tree object
#'
#' @export
#'
setMethod(
	'subtree',
	signature(
		x = 'lineage_tree'
	),
	function(
		x,
		leaves = NULL,
		...
	){

		stopifnot(!is.null(leaves))

		is_leaf <- degree(x@graph, mode = 'out') == 0
		stopifnot(all(leaves %in% names(x@x)[is_leaf]))
	
		d <- distances(x@graph, leaves, mode = 'in')
	  v <- names(x@x)[(!is.infinite(d)) %>% colSums() > 0]  # subgraphs that connect to the leaves
		g <- induced_subgraph(x@graph, v)
		x@x <- x@x[v]
		x@graph <- g
		x
	}
) # subtree


#' subtree
#'
#' Extract a subtree with specific leaves
#'
#' @param x a phylo object
#' @param leaves leaves of the extracted tree
#' @param ... additional parameters
#'
#' @return a pylo object
#'
#' @export
#'
setMethod(
	'subtree',
	signature(
		x = 'phylo'
	),
	function(
		x,
		leaves = NULL,
		...
	){
		x <- x %>% as_igraph()
		d <- distances(x, leaves, mode = 'in')
	  v <- as.character(V(x)$name)[(!is.infinite(d)) %>% colSums() > 0]  # subgraphs that connect to the leaves
		x <- induced_subgraph(x, v)
		x %>% prune() %>% as_phylo()
	}
) # subtree


