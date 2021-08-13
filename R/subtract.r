#' subtract
#'
#' Subtract a subtree from a large tree
#'
#' @param x a lineage_tree object
#' @param y a lineage_tree object
#' @param ... additional parameters
#'
#' @return a lineage_tree object
#'
#' @export
#'
setMethod(
	'subtract',
	signature(
		x = 'lineage_tree',
		y = 'lineage_tree'
	),
	function(
		x,
		y,
		...
	){

		# need to verify if y is a subtree of x

		is_leaf <- degree(x@graph, mode = 'out') == 0
		leaves_x <- names(x@x[is_leaf])

		is_leaf <- degree(y@graph, mode = 'out') == 0
		leaves_y <- names(y@x[is_leaf])

		stopifnot(all(leaves_y %in% leaves_x))

		leaves <- leaves_x[!leaves_x %in% leaves_y]

		stopifnot(length(leaves) > 0)
		subtree(x, leaves)
	}
) # subtree

