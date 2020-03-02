#' function for sampling k-mers in the simulated data
#'
sample_kmer <- function(
  n_nodes = 10L,  # number of randomly sampled nodes (leaves or internval nodes)
  k = 1L,  # k-mer,
  alphabets = NULL,
  division = 20L,
  sequence_length = 200L,  # sequence length
  mutation_prob = 0.1,  # mutation proability
  outcome_prob = NULL,
  deletion = TRUE  # whether or not including the deletion events
){

  states <- do.call('paste0', do.call('expand.grid', lapply(1:k, function(j) alphabets)))                                                          
  sim <- simulate(
    n_samples = 200L,     # number of samples to simulate
    n_targets  = sequence_length,  # number of targets
    mutation_prob = mutation_prob,   # mutation proability
    division = division,        # number of cell divisons
    alphabets = alphabets, # possible mutational outcomes
    outcome_prob = outcome_prob,  # outcome probability of each letter
    deletion = deletion           # whether or not considering inter-mutatoin deletion
  )
  s <- sim@x %>%
    length() %>%
    sample(n_nodes)
  X <- sim@x[s] %>% as.character()
  p <- expand.grid(
    from = 1:n_nodes, 
    to = 1:n_nodes, 
    start = 1:(sim@n_targets - k + 1)
    ) %>%
    filter(from < to)
  D <- sim@graph %>% distances(
    v = names(sim@x)[s], 
    to = names(sim@x)[s]
  )
  d <- D[cbind(p[, 'from'], p[, 'to'])]
  # the k-mer at each position for each pairs of sample
  str_from <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'from'], p[, 'start'] + j - 1)]))
  str_to <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'to'], p[, 'start'] + j - 1)]))
  data.frame(
    from = factor(str_from, states), 
    to = factor(str_to, states), 
    distance = d
  )
}   
