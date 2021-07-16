#' summarize_kmer
#'
#' Summarize kmer distributions with input sequences
#'
#' @param x input data as a phyDat object
#' @param division number of cell division
#' @param k k-mer (default = 2)
#' @param reps number of simulated trees
#' @param n_samples number of samples to simulate
#' @param n_nodes number of nodes to sample (including both leaves and internval nodes)
#' @param n_targets number of targets
#' @param n_targets sequence length. If this argument is missing, the length of the input sequences will be used.
#'
#' @return a kmer_summary object
#'
#' @author Wuming Gong (gongx030@umn.edu)
#' 
#' @export
#'
setMethod(
  'summarize_kmer',
  signature(
    x = 'phyDat'
  ),
  function(
    x,
    division = 16L,
    k = 2,
    reps = 20L,
		n_samples = 200L,
		n_nodes = 100L,
		n_targets 
  ){

		config <- process_sequence(x, division = division)

		if (!missing(n_targets))
			config@n_targets <- n_targets

    summarize_kmer_core(k, reps, n_samples, n_nodes, config)
  }
)


#' summarize_kmer_core
#'
#' Summarize kmer distributions (core function)
#'
#' @param k k-mer (default = 2)
#' @param reps number of simulated trees
#' @param n_samples number of samples to simulate
#' @param n_nodes number of nodes to sample (including both leaves and internval nodes)
#' @param config lineage tree configuration (a lineage_tree_config object)
#'
#' @return a kmer_summary object
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
summarize_kmer_core <- function(
  k = 2,
  reps = 20L,
	n_samples = 200L,
	n_nodes = 100L,
	config = NULL
){

	alphabets <- config@alphabets[config@frequency > 0]
	max_distance <- (config@division - 1) * 2
	kmers <- do.call('paste0', do.call('expand.grid', lapply(1:k, function(j) alphabets)))

  sprintf('simulating | k=%d | alphabets=%d | mutation=%.3f | division=%d | sample=%d | n_targets=%d | dropout_prob=%.3f', k, length(alphabets), config@mutation_prob, config@division, reps, config@n_targets, config@dropout_prob) %>%
		message()

  p <- expand.grid(
    from = 1:n_nodes,
    to = 1:n_nodes,
    start = 1:(config@n_targets - k + 1)
   ) %>%
    filter(.data$from < .data$to)

#  df <- do.call('rbind', bplapply(
  df <- do.call('rbind', lapply(
    seq_len(reps), 
    function(i){
      sim <- simulate(n_samples = n_samples, config)    
      s <- sample(length(sim@x), n_nodes)
      X <- sim@x[s] %>% as.character()
      D <- sim@graph %>% distances(
        v = names(sim@x)[s],
        to = names(sim@x)[s]
      )
      d <- D[cbind(p[, 'from'], p[, 'to'])]

      str_from <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'from'], p[, 'start'] + j - 1)]))
      str_to <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'to'], p[, 'start'] + j - 1)]))

      data.frame(
        from = str_from,
        to = str_to,
        distance = d
      ) %>%
        group_by(.data$from, .data$to, .data$distance) %>%
        tally()
    }
  )) %>%
    group_by(.data$from, .data$to, .data$distance) %>%
    summarize(n = sum(n))

  new(
    'kmer_summary',
    df = df,
		k = k,
    reps = reps,
		alphabets = alphabets,
		max_distance = max_distance,
		kmers = kmers,
		config = config
  )
} # summarize_kmer_core


