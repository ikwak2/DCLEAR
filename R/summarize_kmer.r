setGeneric('summarize_kmer', function(x, ...) standardGeneric('summarize_kmer'))

#' summarize_kmer
#'
#' Summarize kmer distributions with input sequences
#' @param x input data as a phyDat object
#'
#' @return a kmer_summary object
#'
#' @author Wuming Gong (gongx030@umn.edu)
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
    mutation_prob
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

    if (missing(mutation_prob))
      mutation_prob <- 1 - (freq[DEFAULT] / sum(freq))^(1 / division)

    freq[DEFAULT] <- 0
    outcome_prob <- freq / sum(freq)

    alphabets <- names(outcome_prob)

    sequence_length <- ncol(x %>% as.character())

    summarize_kmer_core(k, alphabets, mutation_prob, division, reps, sequence_length, outcome_prob)

  }
)


#' summarize_kmer
#'
#' Summarize kmer distributions without input sequences
#'
#' @return a kmer_summary object
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
  'summarize_kmer',
  signature(
    x = 'missing'
  ),
  function(
    x,
    division = 16L,
    k = 2,
    reps = 20L,
    mutation_prob = 0.1,
    alphabets,
    sequence_length = 200L
  ){
    DEFAULT <- '0'
    DELETION <- '-'
    DROPOUT <- '*'

    outcome_prob <- rowMeans(do.call('cbind', bplapply(seq_len(1000L), function(i) sample_outcome_prob(alphabets))))

    summarize_kmer_core(k, alphabets, mutation_prob, division, reps, sequence_length, outcome_prob)
  }
)

#' summarize_kmer_core
#'
#' Summarize kmer distributions (core function)
#'
#' @param k k-mer (default = 2)
#' @param alphabets alphabets used in the tree
#' @param mutation_prob mutation probability
#' @param division number of cell division
#' @param reps number of simulated trees
#' @param sequence_length sequence length (e.g. number of targets)
#' @param outcome_prob outcome probability of each letter
#' @param n_nodes number of sampled tip or internval nodes in the simulated tree (default: 100L)
#'
#' @return a kmer_summary object
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
summarize_kmer_core <- function(
  k = 2,
  alphabets,
  mutation_prob,
  division,
  reps,
  sequence_length,
  outcome_prob,
  n_nodes = 100L
){

  kmers <- do.call('paste0', do.call('expand.grid', lapply(1:2, function(j) alphabets)))

  max_distance <- (division - 1) * 2

  flog.info(sprintf('simulating | k=%d | alphabets=%d | mutation=%.3f | division=%d | sample=%d | sequence_length=%d', k, length(alphabets), mutation_prob, division, reps, sequence_length))

  p <- expand.grid(
    from = 1:n_nodes,
    to = 1:n_nodes,
    start = 1:(sequence_length - k + 1)
   ) %>%
    filter(from < to)

  df <- do.call('rbind', bplapply(
    seq_len(reps), 
    function(i){

      sim <- simulate(
        n_samples = 200L,     # number of samples to simulate
        n_targets  = sequence_length,  # number of targets
        mutation_prob = as.numeric(mutation_prob),   # mutation proability
        division = division,        # number of cell divisons
        alphabets = alphabets, # possible mutational outcomes
        outcome_prob = as.numeric(outcome_prob),  # outcome probability of each letter
        deletion = TRUE # whether or not considering inter-mutatoin deletion
      ) 

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
        group_by(from, to, distance) %>%
        tally()
    }
  )) %>%
    group_by(from, to, distance) %>%
    summarize(n = sum(n))

  flog.info('finished')


  new(
    'kmer_summary',
    df = df,
    alphabets = alphabets,
    kmers = kmers,
    outcome_prob = as.numeric(outcome_prob),
    sequence_length = sequence_length,
    division = division,
    mutation_prob = as.numeric(mutation_prob),
    reps = reps,
    max_distance = 2 * (division - 1),
    k = k
  )
} # summarize_kmer_core
