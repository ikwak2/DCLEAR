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

		tree <- random_tree(n_samples = n_samples, division = config@division) %>%
			as.data.frame()

		mutation_site <- sample_mutation_site(tree, config)
		outcome <- sample_mutation_outcome(tree, mp = NULL, config)

		simulate_core(config, tree, mutation_site, outcome)

	}
) # simulate


#' simulate
#'
#' Simulate a cell lineage tree based on a set of sequences
#'
#' @param config simulation configuration; a lineage_tree_config object
#' @param x a sequence object
#' @param n_samples number of samples to simulate
#' @param k Number of trials
#' @param greedy Whether ot not use a greedy search
#' @param ... additional parameters
#'
#' @return a lineage_tree object
#'
#' @author Wuming Gong (gongx030@umn.edu)
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
		k = 50,
		greedy = TRUE,
		...
	){

		# compute site specific outcome probability
		mp <- x %>% as.character() %>% positional_mutation_prob(config)
		tree <- random_tree(n_samples = n_samples, division = config@division) %>%
			as.data.frame()

		if (greedy){

			for (h in 1:config@division){
				tree_h <- tree[tree$height <= h, ]
				ms <- lapply(1:k, function(i) sample_mutation_site(tree_h, config))
				oc <- lapply(1:k, function(i) sample_mutation_outcome(tree_h, mp, config))
				for (i in 1:k){
					if (h > 1){
						ms[[i]][rownames(mutation_site), ] <- mutation_site
						oc[[i]][rownames(outcome), ] <- outcome
					}
				}
				raw <- lapply(1:k, function(i) get_sequence(ms[[i]], tree_h, oc[[i]], config))
				scores <- sapply(1:k, function(i) x %>% as.character() %>% score_simulation(raw[[i]], config))
				sprintf('h=%d | best score=%.3f', h, max(scores)) %>% message()
				i <- which.max(scores)
				mutation_site <- ms[[i]]
				outcome <- oc[[i]]
			}
		}else{
			mutation_site <- sample_mutation_site(tree, config)
			outcome <- sample_mutation_outcome(tree, mp, config)
			raw <- get_sequence(mutation_site, tree, outcome, config)
		}

		simulate_core(config, tree, mutation_site, outcome)

	}
)


#' simulate_core
#'
#' Simulate a cell lineage tree
#' Adoped from https://github.com/elifesciences-publications/CRISPR_recorders_sims/blob/master/MATLAB_sims/GESTALT_30hr_1x_simulation.m
#'
#' @param config simulation configuration; a lineage_tree_config object
#' @param tree a matrix representing the lineage tree
#' @param mutation_site a binary matrix indicating the mutation sites
#' @param outcome a character matrix
#' @return a `lineage_tree` object
#'
simulate_core <- function(config, tree, mutation_site, outcome){

	raw <- get_sequence(mutation_site, tree, outcome, config)

	if (config@deletion){
		x <- add_deletion(raw, tree, mutation_site, config)
	}else{
		x <- raw
	}

	if (config@dropout_prob > 0){
		x <- add_dropout(x, config)
	}

	x <- x %>% phyDat(type = 'USER', levels = config@alphabets)
	raw <- raw %>% phyDat(type = 'USER', levels = config@alphabets)

	new('lineage_tree', 
		x = x, 
		raw = raw,
		graph = as_igraph(tree, config), 
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
		apply(2, function(y) table(factor(y, config@alphabets)))
  mp[config@outcome_prob == 0, ] <- 0
	w <- 1 / colSums(mp)
	w[is.infinite(w)] <- 0
	mp <- mp %*% diag(w)
	mp[config@default_character, w == 0] <- 1
	mp
}


#' sample_mutation_site
#' 
#' Sample mutation site
#'
#' @param tree a data frame
#' @param config a lineage_tree_config object
#' @return a mutation site matrix
#'
sample_mutation_site <- function(tree, config){

	ancestor <- 1
	ms <- matrix(FALSE, nrow = 1, ncol = config@n_targets, dimnames = list(ancestor %>% get_node_names(), NULL))	# mutation sites
	is_default <- matrix(TRUE, nrow = 1, ncol = config@n_targets, dimnames = list(ancestor %>% get_node_names(), NULL))	
	h <- 1	# a unit time
	while (h < config@division){
		pairs <- tree %>% filter(.data$height == h) 	# parent-child pairs at level h
		parents <- pairs[, 1] %>% get_node_names()
		children <- pairs[, 2] %>% get_node_names()
		d <- is_default[parents, , drop = FALSE]	# whether or not states in parent cells are default
		rownames(d) <- children
		mutation_site <- runif(config@n_targets * nrow(pairs)) < config@mutation_prob & d # randomly sample the sites are mutable and default in parents
		d[mutation_site] <- FALSE	
		ms <- rbind(ms, mutation_site)
		is_default <- rbind(is_default, d)
		h <- h + 1
	}
	ms
}

#' sample_mutation_outcome
#'
#' Sample mutation outcome
#'
#' @param x an igraph object
#' @param mp a mutation site matrix
#' @param config a lineage_tree_config object
#' @return a outcome matrix
#'
sample_mutation_outcome <- function(x, mp = NULL, config){

	x <- as_igraph(x, config)

	if (is.null(mp)){
		mp <- matrix(
			config@outcome_prob, 
			nrow = length(config@outcome_prob), 
			ncol = config@n_targets, 
			dimnames = list(config@alphabets, NULL)
		)
	}

	outcome <- matrix(config@default_character, nrow = length(V(x)), ncol = config@n_targets, dimnames = list(V(x)$name, NULL))
	for (j in 1:config@n_targets){	# for each position
		outcome[, j] <- sample(config@alphabets, nrow(outcome), prob = mp[, j], replace = TRUE)
	}
	outcome
}


#' get_sequence
#'
#' Get sequencees
#'
#' @param x a character matrix
#' @param tree a matrix representing the lineage tree
#' @param outcome a character matrix
#' @param config a lineage_tree_config object
#' @return a character matrix
#'
get_sequence <- function(x, tree, outcome, config){
	ancestor <- 1
	h <- 1	# a unit time
	y<- matrix(config@default_character, nrow = 1L, ncol = config@n_targets, dimnames = list(ancestor %>% get_node_names(), NULL))
	while (h < config@division){
		pairs <- tree %>% filter(.data$height == h) 	# parent-child pairs at level h
		parents <- pairs[, 1] %>% get_node_names()
		children <- pairs[, 2] %>% get_node_names()
		r <- y[parents, ]
		rownames(r) <- children
		m <- x[children, ]
		r[m] <- outcome[children, ][m]
		y <- rbind(y, r)
		h <- h + 1
	}
	y
}


#' add_deletion
#' 
#' Add deletion
#' @param x a character matrix
#' @param tree a matrix representing the lineage tree
#' @param mutation_site a binary matrix for mutation site
#' @param config a lineage_tree_config object
#' @return a character matrix with deletions
#'
add_deletion <- function(x, tree, mutation_site, config){

	h <- 1	# a unit time
	while (h < config@division){
		pairs <- tree %>% filter(.data$height == h) 	# parent-child pairs at level h
		parents <- pairs[, 1] %>% get_node_names()
		children <- pairs[, 2] %>% get_node_names()

		for (i in 1:nrow(pairs)){ # for each child
			del <- x[parents[i], ] == config@deletion_character	# pass the deletion marks on parents
			x[children[i], del] <- config@deletion_character
			ms <- mutation_site[children[i], ] & x[children[i], ] != config@deletion_character
			n_mut <- sum(ms)
			if (n_mut > 1){
				# find all pairs of mutation sites within 20 bp
				p <- expand.grid(from = which(ms), to = which(ms)) %>%
					filter(to - from <= 20 & to - from >= 2)
				if (nrow(p) > 0){
					j <- sample(1:nrow(p), 1)
					from <- p[j, 'from'] + 1
					to <- p[j, 'to'] - 1
					x[children[i], from:to] <- config@deletion_character
				}
			}
		}
		h <- h + 1
	}
	x
}

#' add_dropout
#' 
#' Add dropout events
#' 
#' @param x a character matrix
#' @param config a lineage_tree_config object
#' @return a character matrix with dropout events
#' 
add_dropout <- function(x, config){
	# randomly add dropout events
	dropout_position <- sample(prod(dim(x)), round(prod(dim(x)) * config@dropout_prob), replace = FALSE)
	x[dropout_position] <- config@dropout_character
	x
}

#' score_simulation 
#' 
#' Compare two sets of sequences
#'
#' @param x a character matrix
#' @param y a character matrix
#' @param config a lineage_tree_config object
#' @return numeric scores
#'
score_simulation <- function(x, y, config){
	p <- expand.grid(from = 1:nrow(x), to = 1:nrow(y))
	s <- matrix(rowSums(x[p[, 'from'], ] == y[p[, 'to'], ]), nrow = nrow(x), ncol = nrow(y))
	mean(rowMaxs(s) / rowSums(x != config@dropout_character & x != config@deletion_character))
}

