setGeneric('compute_replacement_matrix', function(x, ...) standardGeneric('compute_replacement_matrix'))

#' compute_replacement_matrix
#'
#' Simulate replacement matrix with input sequence data
#'
setMethod(
	'compute_replacement_matrix',
	signature(
		x = 'phyDat'
	),
	function(
		x,
		n_batch = 10L,
		n_cells_per_batch = 200L,
		k = 2,
		mc.cores = 2,
		division = 20L,
    dropout = TRUE
	){

		DEFAULT <- '0'
	  DELETION <- '-'
	  DROPOUT <- '*'

		alphabets <- levels(x)
		freq <- x %>% 
			as.character() %>% 
			factor(alphabets) %>% 
			table()
		freq <- freq[freq > 0]

		if (any(names(freq) %in% DELETION))
			freq[DELETION] <- 0

		if (any(names(freq) %in% DROPOUT))
			freq[DROPOUT] <- 0

		mutation_probs <- 1 - (freq[DEFAULT] / sum(freq))^(1 / division)
		freq[DEFAULT] <- 0
		outcome_prob <- freq / sum(freq)
		alphabets <- names(outcome_prob)
		sequence_length <- ncol(x %>% as.character())

		compute_replacement_matrix_core(
			n_batch = n_batch,
			n_cells_per_batch = n_cells_per_batch,
			k = k,
			alphabets = alphabets,
			outcome_prob = as.numeric(outcome_prob),
			sequence_length = sequence_length,
			mc.cores = mc.cores,
			mutation_probs = mutation_probs,
			division = division,
      dropout = dropout
		)
	}
) # compute_replacement_matrix



#' compute_replacement_matrix_core
#' 
#' @param n_batch number of simulated trees
#' @param n_cells_per_batch number of sampled tip or internval nodes in the simulated tree
#' @param k k-mer
#' @param alphabets alphabets used in the tree
#' @param outcome_prob outcome probability of each letter
#' @param sequence_length sequence length (e.g. number of targets)
#'
compute_replacement_matrix_core <- function(
	n_batch = 10L,
	n_cells_per_batch = 200L,
	k = 2,
	alphabets = NULL,
	outcome_prob = NULL,
	sequence_length = 200L,
	mutation_probs = seq(0.01, 0.3, by = 0.01),
	division = 20L,
	n_leaves = 200L,
  dropout = TRUE,
	mc.cores = 2
){

	res <- mclapply(1:n_batch, function(i){

		# sampling a mutation proability
		mp <- sample(mutation_probs, 1)	

		flog.info(sprintf('simulating | sample=%4.d/%4.d | k=%2.d | mutation prob=%.3f', i, n_batch, k, mp))

		sim <- simulate(
			n_samples = n_leaves, 		# number of samples to simulate
			n_targets  = sequence_length,  # number of targets
			mutation_prob = mp,		# mutation proability
 		  division = division, 				# number of cell divisons
		  alphabets = alphabets, # possible mutational outcomes
			outcome_prob = outcome_prob,	# outcome probability of each letter
			dropout = dropout 					# whether or not considering inter-mutatoin deletion
		)
	
		# randomly sample nodes, including either leaves or internal nodes
		s <- sample(1:length(sim@x), n_cells_per_batch)	
		X <- sim@x[s] %>% as.character()
		p <- expand.grid(
			from = seq_len(n_cells_per_batch),
			to = seq_len(n_cells_per_batch),
			start = 1:(sequence_length - k + 1)
		) %>% 
			filter(from < to)

		# the shortest distance matrix between any two sampled nodes
		D <- sim@graph %>% 
			distances(v = names(sim@x)[s], to = names(sim@x)[s])

		# the k-mer at each position for each pairs of sample
		str_from <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'from'], p[, 'start'] + j - 1)]))
		str_to <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'to'], p[, 'start'] + j - 1)]))

		# compute the mean distances between any pairs of k-mers
		d <- D[cbind(p[, 'from'], p[, 'to'])]
		d_mean <- aggregate(d, list(from = str_from, to = str_to), mean)

		d_mean

	}, mc.cores = mc.cores)

	# compute the mean distance between any pairs of k-mers over all simulated samples
	r <- do.call('rbind', lapply(1:n_batch, function(i) res[[i]]))
	r <- aggregate(r$x, list(from = r$from, to = r$to), mean)

	# all possible states
	states <- do.call('paste0', do.call('expand.grid', lapply(1:k, function(j) alphabets)))

	# A replacment matrix between two types of k-mers
	S <- sparseMatrix(
		i = as.numeric(factor(r$from, states)), 
		j = as.numeric(factor(r$to, states)), 
		x = r$x, dims = c(length(states), length(states)), 
		dimnames = list(states, states)
	)
	S <- S + t(S)
	S
} # compute_replacement_matrix_core



