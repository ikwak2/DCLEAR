#' simulate
#'
#' Simulate a cell lineage tree
#' Adoped from https://github.com/elifesciences-publications/CRISPR_recorders_sims/blob/master/MATLAB_sims/GESTALT_30hr_1x_simulation.m
#'
#' @param config simulation configuration; a lineage_tree_config object
#' @param x missing 
#' @param n_samples number of samples to simulate
#' @param ... additional parameters
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
		config = 'lineage_tree_config',
		x = 'missing'
	),
	function(
		config,
		x,
		n_samples = 200, # number of samples to simulate 
		...
	){

		simulate_core(config, mp = NULL, n_samples = n_samples, ...)

	}
) # simulate


#' simulate
#'
#' Simulate a cell lineage tree based on a set of sequences
#'
#' @param config simulation configuration; a lineage_tree_config object
#' @param x a sequence object
#' @param n_samples number of samples to simulate
#' @param ... additional parameters
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
		config = 'lineage_tree_config',
		x = 'phyDat'
	),
	function(
		config,
		x,
		n_samples = 200L,
		...
	){

		# compute site specific outcome probability
		mp <- positional_mutation_prob(x, config)
		simulate_core(config, mp = mp, n_samples = n_samples, ...)
	}
)


#' simulate_core
#'
#' Simulate a cell lineage tree
#' Adoped from https://github.com/elifesciences-publications/CRISPR_recorders_sims/blob/master/MATLAB_sims/GESTALT_30hr_1x_simulation.m
#'
#' @param config simulation configuration; a lineage_tree_config object
#' @param mp site specific mutation probability
#' @param n_samples number of samples to simulate
#' @param ... additional parameters
#'
simulate_core <- function(config, mp = NULL, n_samples = 200L, ...){

	tree <- random_tree(n_samples = n_samples, division = config@division) %>%
		as.data.frame()

	if (is.null(mp)){
		mp <- matrix(
			config@outcome_prob, 
			nrow = length(config@outcome_prob), 
			ncol = config@n_targets, 
			dimnames = list(config@alphabets, NULL)
		)
	}

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
				xc[i, mutation_site] <- sapply(which(mutation_site), function(j) sample(config@alphabets, 1L, prob = mp[, j]))
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

	leaf_min <- 2^(division - 1)
	leaf_max <- (2^division - 1)
	to <- (sample.int(leaf_max - leaf_min + 1L, n_samples) + leaf_min - 1L) %>% 
		sort()
	h <- division - 1	# leaf level
	edges <- NULL
	while (h > 0){
		from <- floor(to / 2)
		edges <- rbind(cbind(from = from, to = to, height = h), edges)
		to <- sort(unique(from))
		h <- h - 1
	}
	edges
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
get_node_names <- function(x) sprintf('node_%s', sprintf('%.0f', x) %>% str_pad(20, pad = '0'))


#' positional_mutation_prob
#'
#' Convenient function for get node names
#'
#' @param x a phyDat object
#' @param config a lineage_tree_config object
#'
#' @return a positional mutation probability matrix
#'
#' @export
#'
positional_mutation_prob <- function(x, config){
	mp <- x %>% 
		as.character() %>% 
		apply(2, function(y) table(factor(y, config@alphabets)))
  mp[config@outcome_prob == 0, ] <- 0
	w <- 1 / colSums(mp)
	w[is.infinite(w)] <- 0
	mp <- mp %*% diag(w)
	mp[config@default_character, w == 0] <- 1
	mp
}

