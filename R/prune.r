#' prune
#'
#' Trim a full lineage tree into phylogenetic tree
#'
#' @param x a lineage_tree object
#' @param ... additional parameters passed to as_phylo()
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

		x@graph <- x@graph %>% prune(...)
		x@x <- x@x[V(x@graph)$name %>% unlist()]
		x

	}
)

#' prune
#'
#' Trim a full lineage tree into phylogenetic tree
#'
#' @param x an igraph object 
#' @param weighted whether or not keep the edge weight (default: TRUE)
#' @param ... additional parameters 
#'
#' @return an igraph object 
#'
#' @export
#'
setMethod(
	'prune',
	signature(
		x = 'igraph'
	),
	function(
		x,
		weighted = TRUE,
		...
	){

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
		if (weighted){
			edges <- (x[] %>% summary())[, 1:2] %>% as.matrix()
			E(x)$weight <- d_lb[edges]
		}
		x

	}
)
