setGeneric('dist_replacement', function(x, kmer_summary, ...) standardGeneric('dist_replacement'))

setMethod(
	'dist_replacement',
	signature(
		x = 'phyDat',
    kmer_summary = 'kmer_summary'
	),
	function(x, kmer_summary, method, ...){

    if (method == 'sampling'){
  		dist_replacement_core(x, kmer_summary, ...)
    }else if (method == 'bayesian'){
      dist_replacement_bayesian_core(x, kmer_summary, ...)
    }
  }
)


#' dist_replacement_core
#' 
#' Compute the replacement distance matrix using a replacement matrix derived by sampling
#'
dist_replacement_core <- function(x, kmer_summary, ...){

  p_D_AB <- get_distance_probability(kmer_summary)

  max_distance <- kmer_summary@max_distance
  kmers <- kmer_summary@kmers

  sequence_length <- ncol(x %>% as.character())

	p <- expand.grid(
		from = 1:length(x),
		to = 1:length(x),
    distance = 1:max_distance
	) %>% 
		filter(from < to)

  # input sequence in BStringSet object
  y <- do.call('paste0', as.data.frame(x %>% as.character())) %>%
    BStringSet()

  log_prob <- rep(0, length(x) * (length(x) - 1) / 2 * max_distance)

  for (start in 1:sequence_length){ # for each position in the target sequence

    if (start %% 10 == 0)
      flog.info(sprintf('posterior probability | position=%5.d/%5.d', start, sequence_length))

 		yi <- substr(y, start, start)

    str_from <- yi[p[, 'from']] %>% 
      factor(kmers) %>% 
      as.numeric()

    str_to <- yi[p[, 'to']] %>% 
      factor(kmers) %>% 
      as.numeric()

    log_prob_i <- log(p_D_AB[cbind(str_from, str_to, p$distance)] + 1e-1)
    log_prob <- log_prob + log_prob_i
  }

  P <- matrix(log_prob, nrow = length(x) * (length(x) - 1) / 2, ncol = max_distance)
  post <- exp(P - apply(P, 1, logSumExp))

  p1 <- p %>% filter(distance == 1)
  D <- sparseMatrix(
    i = p1[, 'from'], 
    j = p1[, 'to'], 
    x = rowSums(post %*% Diagonal(x = 1:max_distance)),
    dims = c(length(x), length(x)), 
    dimnames = list(names(x), names(x))
  )
  D <- t(D) + D
  D %>% as.dist()

} # dist_replacement_core

 
#' dist_replacement_bayesian_core
#'
dist_replacement_bayesian_core <- function(x, kmer_summary){

  # prior distribution of nodal distance
  p_D <- kmer_summary  %>% get_distance_prior()

  p_B_AD <- kmer_summary %>% get_conditional_replacement_probability()

  max_distance <- length(p_D)
  alphabets <- dimnames(p_B_AD)[[1]]

  sequence_length <- ncol(x %>% as.character())

	p <- expand.grid(
		from = 1:length(x),
		to = 1:length(x),
    distance = 1:max_distance
	) %>% 
		filter(from < to)

  # input sequence in BStringSet object
  y <- do.call('paste0', as.data.frame(x %>% as.character())) %>%
    BStringSet()

  log_prob <- rep(0, length(x) * (length(x) - 1) / 2 * max_distance)

  for (start in 1:sequence_length){ # for each position in the target sequence

    if (start %% 10 == 0)
      flog.info(sprintf('posterior probability | position=%5.d/%5.d', start, sequence_length))

 		yi <- substr(y, start, start)

    str_from <- yi[p[, 'from']] %>% 
      factor(alphabets) %>% 
      as.numeric()

    str_to <- yi[p[, 'to']] %>% 
      factor(alphabets) %>% 
      as.numeric()

    log_prob_i <- log(p_B_AD[cbind(str_from, str_to, p$distance)] + 1e-10)
    log_prob <- log_prob + log_prob_i
  }

  log_prob <- log_prob + log(p_D[p$distance])

  P <- matrix(log_prob, nrow = length(x) * (length(x) - 1) / 2, ncol = max_distance)
  post <- exp(P - apply(P, 1, logSumExp))

  p1 <- p %>% filter(distance == 1)
  D <- sparseMatrix(
    i = p1[, 'from'], 
    j = p1[, 'to'], 
    x = rowSums(post %*% Diagonal(x = 1:max_distance)),
    dims = c(length(x), length(x)), 
    dimnames = list(names(x), names(x))
  )
  D <- t(D) + D
  D %>% as.dist()

} # dist_replacement_bayesian_core


#' get_distance_prior 
#'
#' prior distribution of distance
#'
get_distance_prior <- function(x){

  d <- x@df %>%
    ungroup() %>%
    group_by(distance) %>%
    summarize(n = sum(n)) %>%
    mutate(prob = n / sum(n)) %>%
    select(distance, prob)

  prior <- rep(0, x@max_distance)
  prior[d$distance] <- d$prob
  prior

} # get_distance_prior


get_conditional_transition_probability <- function(x){

  # conditional transition probability
  d <- x@df %>%
    ungroup() %>%
    filter(k == 2) %>%
    mutate(
      from2 = paste0(substr(from, 1, 1), substr(to, 1, 1)),
      to2 = paste0(substr(from, 2, 2), substr(to, 2, 2))
    ) %>%
    mutate(from = from2, to = to2) %>%
    select(from, to, distance, n) %>%
    right_join(
      expand.grid(
        from = x@dimers, 
        to = x@dimers, 
        distance = 1:x@max_distance, 
        stringsAsFactors = FALSE
      ),
      by = c('from', 'to', 'distance')
    ) %>%
    replace_na(list(n = 10))  %>%
    group_by(distance, from) %>%
    mutate(prob = n / sum(n)) %>%
    select(from, to, distance, prob)

  array(
    d$prob, 
    dim = c(length(x@dimers), length(x@dimers), x@max_distance),
    dimnames = list(x@dimers, x@dimers, NULL)
  )

} # get_conditional_transition_probability


#' get_conditional_replacement_probability
#'
get_conditional_replacement_probability <- function(x){

  # conditional probability of single nucleotide replacement
  d <- x@df %>%
    ungroup() %>%
    select(from, to, distance, n) %>%
    right_join(
      expand.grid(
        from = x@alphabets, 
        to = x@alphabets, 
        distance = 1:x@max_distance, 
        stringsAsFactors = FALSE
      ),
      by = c('from', 'to', 'distance')
    ) %>%
    replace_na(list(n = 0))  %>%
    mutate(n = ifelse(from == to, n + 1, n)) %>%
    group_by(from, distance) %>%
    mutate(prob = n / sum(n)) %>%
    select(from, to, distance, prob)

  array(
    d$prob, 
    dim = c(length(x@alphabets), length(x@alphabets), x@max_distance),
    dimnames = list(x@alphabets, x@alphabets, NULL)
  )

} # get_conditional_replacement_probability

get_letter_prior <- function(x){

  d <- x@df %>%
    ungroup() %>%
    filter(k == 1) %>%
    select(from, n) %>%
    right_join(
      data.frame(
        from = x@alphabets
      ),
      by = 'from'
    ) %>%
    group_by(from) %>%
    summarize(n = sum(n)) %>%
    mutate(prob = n / sum(n))
  d$prob

} # get_letter_prior


#' get_replacement_probability
#'
get_replacement_probability <- function(x){

  d <- x@df %>%
    ungroup() %>%
    filter(k == 1) %>%
    select(from, to, n) %>%
    right_join(
      expand.grid(
        from = x@alphabets, 
        to = x@alphabets, 
        stringsAsFactors = FALSE
      ),
      by = c('from', 'to')
    ) %>%
    replace_na(list(n = 0))  %>%
    group_by(from, to) %>%
    summarize(n = sum(n)) %>%
    mutate(prob  = n / sum(n))

  matrix(
    d$prob, 
    nrow = length(x@alphabets), 
    ncol = length(x@alphabets),
    dimnames = list(x@alphabets, x@alphabets)
  )

} # get_replacement_probability


#' get_distance_probability
#'
#' Compute the conditional distribution of distance according to k-mer replacement
#'
get_distance_probability <- function(x){

  d <- x@df %>%
    ungroup() %>%
    select(from, to, distance, n) %>%
    mutate(n = ifelse(from == to, n + 1, n)) %>%
    right_join(
      expand.grid(
        from = x@kmers, 
        to = x@kmers, 
        distance = 1:x@max_distance, 
        stringsAsFactors = FALSE
      ),
      by = c('from', 'to', 'distance')
    ) %>%
    replace_na(list(n = 0))  %>%
    group_by(from, to) %>%
    mutate(prob  = n / sum(n))

  array(
    d$prob, 
    dim = c(length(x@kmers), length(x@kmers), x@max_distance),
    dimnames = list(x@kmers, x@kmers, NULL)
  )

} # get_distance_probability


