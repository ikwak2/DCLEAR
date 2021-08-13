#' as_igraph
#' 
#' Convert an phylo object to an igraph object, while keeping the weight (in contrast to igraph::as.igraph)
#' @param x a phylo object
#' @return an igraph object
#' 
#' @export
#'
setMethod(
	'as_igraph',
	signature(
		x = 'phylo'
	),
	function(x){
		n <- length(x$tip.label) + x$Nnode
		g <- graph_from_edgelist(x$edge)
		V(g)$name <- c(x$tip.label, sprintf('Node%d', 1:x$Nnode))
		if (!is.null(x$edge.length)){
			E(g)$weight <- x$edge.length
		}
		g
	}
)
