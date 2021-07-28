#' prune
#'
#' Trim a full lineage tree into phylogenetic tree
#'
#' @param x a lineage_tree object
#' @param ... additional parameters
#'
#' @return a list of phyDat and phylo objects
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

		tree <- x@graph %>% as_phylo()
		sequence <- get_leaves(x)
		sequence <- sequence[tree$tip.label]
		list(sequence = sequence, tree = tree)
	}
)
