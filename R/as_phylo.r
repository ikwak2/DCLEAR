#' as_phylo
#' 
#' Convert an igraph object to a phylo object
#' @param x an igraph object
#' @return a phylo object or a igraph object
#' 
#' @export
#'
setMethod(
	'as_phylo',
	signature(
		x = 'igraph'
	),
	function(x){

		# check if this is a phylogenetic tree

		is_leaf <- degree(x, mode = 'out') == 0
		n_leaves <- sum(is_leaf)
		v <- topo_sort(x, mode = 'in')$name %>% 
			rev() %>%
			unlist() 
		v <- c(v[(vcount(x) - n_leaves + 1):vcount(x)], v[1:(vcount(x) - n_leaves)])
		v0 <- V(x)$name %>% as.character()
		x <- permute.vertices(x, factor(v0, v) %>% as.numeric())
		is_leaf <- degree(x, mode = 'out') == 0
		d <- distances(x)
		y <- list()
		y$edge <- (x[] %>% summary())[, 1:2] %>% as.matrix()
		dimnames(y$edge) <- NULL
		y$edge.length <- d[y$edge]
		y$Nnode <- sum(!is_leaf)
		y$tip.label <- V(x)$name[is_leaf] %>% as.character()
		attr(y, 'class') <- 'phylo'
		attr(y, "order") <- 'cladewise'
		y
	}

) # as_phylo


