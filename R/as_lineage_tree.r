#' as_lineage_tree
#' 
#' Convert a phylo object and a phyDat object to a lineage_tree object
#'
#' @param x a phyDat object
#' @param y a phylo object
#' @param config a lineage_tree_config object
#' @param ... additional parameters
#' @return a lineage_tree object
#' 
#' @export
#'
setMethod(
	'as_lineage_tree',
	signature(
		x = 'phyDat',
		y = 'phylo',
		config = 'lineage_tree_config'
	),
	function(x, y, config, ...){
		s <- matrix(config@default_character, y$Nnode, config@n_targets, dimnames = list(sprintf('Node%d', 1:y$Nnode)))
		s <- phyDat(s, type = 'USER', levels = config@alphabets)
		x <- rbind(s, x)
	  y <- y %>% as_igraph()
		x <- x[as.character(V(y)$name)]
		new('lineage_tree',
			x = x,
			graph = y,
			config = config
		)
	}

)
