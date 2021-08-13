#' get_leaves
#'
#' Get the leaf sequences
#'
#' @param x a lineage_tree object
#' @param ... additional parameters
#'
#' @return a phyDat object
#'
#' @export
#'
setMethod(
	'get_leaves',
	signature(
		x = 'lineage_tree'
	),
	function(
		x,
		...
	){
		is_leaf <- degree(x@graph, mode = 'out') == 0
		x@x[is_leaf]
	}
) # get_leaves

