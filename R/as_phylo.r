#' as_phylo
#' 
#' Convert an igraph object to a phylo object
#' @param x an igraph object
#' @return a phylo object
#' 
#' @export
#'
setMethod(
	'as_phylo',
	signature(
		x = 'igraph'
	),
	function(x){

		# reverse the vertic order so that the leaves have the smallest index
		x <- permute.vertices(x, vcount(x):1)	

		is_leaf <- degree(x, mode = 'out') == 0
		is_branch <- degree(x, mode = 'out') == 2
		is_link <- degree(x, mode = 'out') == 1
		d <- distances(x, v = V(x)[is_link], to = V(x)[is_leaf | is_branch], mode = 'out')

		d_lb <- distances(x, v = V(x)[is_leaf | is_branch], to = V(x)[is_leaf | is_branch])

		n_nodes <- sum(is_leaf | is_branch)	# number of nodes in phylo
		new2old <- which(is_leaf | is_branch)	# map from new vertex index to old vertex index
		name2id <- 1:n_nodes
		names(name2id) <- names(new2old)

		old2new <- rep(NA, vcount(x))
		names(old2new) <- V(x)$name
		old2new[is_leaf] <- name2id[names(which(is_leaf))]
		old2new[is_branch] <- name2id[names(which(is_branch))]
		old2new[is_link] <- name2id[names(new2old[apply(d, 1, which.min)])]
		x <- contract(x, old2new)
		x <- simplify(x)
		x <- x %>% set_vertex_attr('name', index = 1:vcount(x), value = names(name2id))
		x <- permute.vertices(x, vcount(x):1)
		is_leaf <- degree(x, mode = 'out') == 0
		y <- list()
		y$edge <- (x[] %>% summary())[, 1:2] %>% as.matrix()
		y$edge.length <- d_lb[y$edge]
		y$Nnode <- sum(is_branch)
		y$tip.label <- V(x)$name[is_leaf]
		attr(y, 'class') <- 'phylo'
		attr(y, "order") <- 'cladewise'
		y
	}

) # as_phylo


