#' as_phylo
#' 
#' Convert an igraph object to a phylo object
#' @param x an igraph object
#' @param weighted whether include the branch length
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
		permutation <- 1:vcount(x)
		permutation[is_leaf] <- which(is_leaf) - n_leaves + 1
		permutation[!is_leaf] <- which(!is_leaf) + n_leaves 
		x <- permute.vertices(x, permutation)
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
		RF.dist(y, y)
		y
	}

) # as_phylo


