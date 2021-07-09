#' simulate
#'
#' Simulate a cell lineage tree
#' Adoped from https://github.com/elifesciences-publications/CRISPR_recorders_sims/blob/master/MATLAB_sims/GESTALT_30hr_1x_simulation.m
#'
#' @param n_samples number of samples to simulate
#' @param config simulation configuration; a lineage_tree_config object
#'
#' @return a lineage_tree object
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
#' @export
#'

setMethod(
	'simulate',
	signature(
		config = 'lineage_tree_config'
	),
	function(
		config,
		n_samples = 200 # number of samples to simulate 
	){

		tree <- random_tree(n_samples = n_samples, division = config@division)

		ancestor <- 1
		x <- matrix(config@default_character, nrow = 1, ncol = config@n_targets, dimnames = list(ancestor %>% get_node_names(), NULL))	# the ancestor sequence
		h <- 1	# a unit time
		while (h < config@division){

			pairs <- tree %>% filter(.data$height == h) 	# parent-child pairs at level h
			parents <- pairs[, 1] %>% get_node_names()
			children <- pairs[, 2] %>% get_node_names()
			xc <- x[parents, , drop = FALSE]	# the states of the children cells
			rownames(xc) <- children

			for (i in 1:nrow(pairs)){ # for each child

				mutation_site <- runif(config@n_targets) < config@mutation_prob & xc[i, ] == config@default_character 	# sampling mutation site
				n_mut <- sum(mutation_site)

				if (n_mut > 0){
					xc[i, mutation_site] <- sample(config@alphabets, n_mut, replace = TRUE, prob = config@outcome_prob)
					if (n_mut > 1 && config@deletion){
						# find all pairs of mutation sites within 20 bp
						p <- expand.grid(from = which(mutation_site), to = which(mutation_site)) %>%
						filter(to - from <= 20 & to - from >= 2)
						if (nrow(p) > 0){
							j <- sample(1:nrow(p), 1)
							from <- p[j, 'from'] + 1
							to <- p[j, 'to'] - 1
							xc[i, from:to] <- config@deletion_character
						}
					}
				}
			}

			x <- rbind(x, xc)
			h <- h + 1
		}

	if (config@dropout_prob > 0){
		# randomly add dropout events
		dropout_position <- sample(prod(dim(x)), round(prod(dim(x)) * config@dropout_prob), replace = FALSE)
		x[dropout_position] <- config@dropout_character
	}

	num_nodes <- 2^config@division - 1
	node_names <- (1:num_nodes) %>% get_node_names()

	A <- sparseMatrix(
		i = tree[, 'from'], 
		j = tree[, 'to'], 
		dims = c(num_nodes, num_nodes),
		dimnames = list(node_names, node_names)
	)
	A <- A[rownames(x), rownames(x)]

	g <- A %>% as('dgCMatrix') %>% graph.adjacency()

	x <- x %>% phyDat(type = 'USER', levels = config@alphabets)

	new('lineage_tree', 
		x = x, 
		graph = g, 
		config = config
	)
	}
) # simulate


#' random_tree
#'
#' Simulate a random lineage tree
#'
#' @param n_samples number of samples to simulate
#' @param division number of cell division
#'
#' @return a data frame
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
random_tree <- function(n_samples, division = 16L){

	ancestor <- 1	# ancestor index
	num_nodes <- 2^division - 1	# number of total nodes in the binary tree
	k <- 2:num_nodes
	edges <- cbind(from = k %/% 2, to = k)	# a binary tree
	edges <- edges[order(edges[, 'to']), ]	# so that edge `from | to` is simply edges[to - 1, ]

	leaves <- 2^(division - 1):(2^division - 1)	# leaf index
	sampled_leaves <- sample(leaves, n_samples)	# leaves to keep

	sampled_edges <- list()
	h <- division - 1	# leaf level
	sampled_edges[[h]] <- edges[sampled_leaves - 1, ]
	while (h > 1){
		to <- unique(sampled_edges[[h]][, 'from'])
		sampled_edges[[h - 1]] <- edges[to - 1, ]
		h <- h - 1
	}

	height <- rep(seq_len(division  - 1), sapply(sampled_edges, nrow))
	sampled_edges <- do.call('rbind', sampled_edges)
	data.frame(sampled_edges, height = height)
} # random_tree


#' sample_outcome_prob
#' 
#' Sampling outcome probability based on a gamma distribution
#'
#' @param config a lineage_tree_config object
#' @param num_states number of states used in simulation.
#' @param shape shape parameter in gamma distribution
#' @param scale scale parameter in gamma distribution
#' 
#' @return a probability vector for each alphabet
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
#' @export
#'
sample_outcome_prob <- function(config, num_states = 20L, shape = 0.1, scale = 2){

	rvs <- rgamma(num_states, shape = shape, scale = scale)
	rvs <- rvs / sum(rvs)
	rvs <- sort(rvs, decreasing = TRUE)

	# the target mutation rate
	config@outcome_prob[4:(4 + num_states - 1)] <- rvs
	names(config@outcome_prob) <- config@alphabets
	config

} # sample_outcome_prob

#' get_node_names
#'
#' Convenient function for get node names
#'
#' @param x node id
#'
#' @return node names
#' @author Wuming Gong (gongx030@umn.edu)
#'
get_node_names <- function(x) sprintf('node_%s', sprintf('%d', x) %>% str_pad(15, pad = '0'))


