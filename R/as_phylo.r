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

		y <- list()	# a phylo object

		# reverse the vertic order so that the leaves have the smallest index
		x <- permute.vertices(x, vcount(x):1)	

		is_leaf <- degree(x, mode = 'out') == 0
		is_branch <- degree(x, mode = 'out') == 2
		is_link <- degree(x, mode = 'out') == 1
		d <- distances(x, v = V(x)[is_link], to = V(x)[is_leaf | is_branch], mode = 'out')

		d_lb <- distances(x, v = V(x)[is_leaf | is_branch], to = V(x)[is_leaf | is_branch])

		n_nodes <- sum(is_leaf | is_branch)
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
			is_leaf <- degree(x, mode = 'out') == 0
		x <- x %>% set_vertex_attr('name', index = 1:vcount(x), value = names(name2id))
		y$edge <- (x[] %>% summary())[, 1:2] %>% as.matrix()
		y$tip.label <- names(name2id)[is_leaf]
		y$Nnode <- vcount(x)
		y$edge.length <- d_lb[y$edge]

		attr(y, 'class') = 'phylo'
		tmpfn <- tempfile(fileext = '.nw')
		write.tree(y, tmpfn)
		read.tree(tmpfn)
	}

) # as_phylo
