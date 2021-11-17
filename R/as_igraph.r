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

#' as_igraph
#' 
#' Convert an phylo object to an igraph object, while keeping the weight (in contrast to igraph::as.igraph)
#' @param x a phylo object
#' @param config a `lineage_tree_config` object
#' @return an igraph object
#' 
#' @export
#'
setMethod(
	'as_igraph',
	signature(
		x = 'data.frame'
	),
	function(x, config){

		v <- unique(c(x[, 'from'], x[, 'to'])) %>% get_node_names()

	  A <- sparseMatrix(
			i = x[, 'from'] %>% get_node_names() %>% factor(v) %>% as.numeric(),
			j = x[, 'to'] %>% get_node_names %>% factor(v) %>% as.numeric(),
			dims = c(length(v), length(v)),
			dimnames = list(v, v),
		) %>% 
			as('dgCMatrix') %>% 
			graph.adjacency()
	}
)
