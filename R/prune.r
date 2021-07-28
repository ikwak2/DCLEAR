#' prune
#'
#' Trim a full lineage tree into phylogenetic tree
#'
#' @param x a lineage_tree object
#' @param ... additional parameters
#'
#' @return a lineage_tree object
#'
#' @export
#'
setMethod(
	'prune',
	signature(
		x = 'lineage_tree'
	),
	function(
		x,
		...
	){

		x@x <- get_leaves(x)
		x@graph %>% as_phylo() %>% as_igraph
		x
	}
)
