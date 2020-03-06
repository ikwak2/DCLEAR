setGeneric('dist_replacement', function(x, replacement_matrix, kmer_summary, ...) standardGeneric('dist_replacement'))

#' dist_replacement
#'
#' compute the kmer replacement distance matrix
#'
#' @export
#
setMethod(
	'dist_replacement',
	signature(
		x = 'phyDat',
		replacement_matrix = 'dgCMatrix',
    kmer_summary = 'missing'
	),
	function(x, replacement_matrix, kmer_summary, ...){
		dist_replacement_core(x, replacement_matrix, ...)
	}
) # 


setMethod(
	'dist_replacement',
	signature(
		x = 'phyDat',
		replacement_matrix = 'missing',
    kmer_summary = 'missing'
	),
	function(x, replacement_matrix, kmer_summary, ...){

 		replacement_matrix <- compute_replacement_matrix(x = x, ...)
 		dist_replacement_core(x, replacement_matrix, ...)
  }
)


setMethod(
	'dist_replacement',
	signature(
		x = 'phyDat',
		replacement_matrix = 'missing',
    kmer_summary = 'kmer_summary'
	),
	function(x, replacement_matrix, kmer_summary, ...){

    prior <- kmer_summary %>% 
      get_distance_prior()

    p_AB_D <- kmer_summary %>%
      get_conditional_replacement_probability()

    p_B_AD <- NULL
    if (kmer_summary@k == 2){
      p_B_AD <- kmer_summary %>% 
       get_conditional_transition_probability()
    }

    dist_replacement_bayesian_core(x, prior, p_AB_D, p_B_AD)

  }
) 


#' dist_replacement_core
#' 
#' Compute the replacement distance matrix using a replacement matrix derived by sampling
#'
dist_replacement_core <- function(x, replacement_matrix){

	X <- x %>% as.character()

	sequence_length <- ncol(X)

	p <- expand.grid(
		from = 1:nrow(X), 
		to = 1:nrow(X)
	) %>% 
		filter(from < to)

	y <- do.call('paste0', as.data.frame(X)) %>%
		BStringSet()

	states <- rownames(replacement_matrix)
	k <- nchar(states)[1]

	D <- matrix(0, nrow = nrow(X), ncol = nrow(X), dimnames = list(names(x), names(x)))	# sum of distance
	res <- NULL

	for (start in 1:c(sequence_length - k + 1)){

		yi <- substr(y, start, start + k - 1)

		str_from <- yi[p[, 'from']] %>% 
			factor(states) %>% 
			as.numeric()

		str_to <- yi[p[, 'to']] %>% 
			factor(states) %>% 
			as.numeric()

		d <- replacement_matrix[cbind(str_from, str_to)]	# distance
																												
		# there are k-mer pairs not present in the simulation
		# * if these are different k-mer's , they have large distance
		# * if these are the same k-mer's, they should be very close
		rare <- d == 0 & str_from != str_to
		d[rare] <- max(d)	# (division - 1) * 2

		Di <- sparseMatrix(
			i = p[, 'from'], 
			j = p[, 'to'], 
			x = d, 
			dims = c(nrow(X), nrow(X)), 
			dimnames = list(names(x), names(x))
			) %>% as.matrix()

		D <- D + Di
	}

	D <- D + t(D)
	diag(D) <- 0
	D %>% as.dist()
}

 
#' dist_replacement_bayesian_core
#'
dist_replacement_bayesian_core <- function(x, prior, p_AB_D, p_B_AD){

  max_distance <- length(prior)
  alphabets <- dimnames(p_AB_D)[[1]]
  dimers <- dimnames(p_B_AD)[[1]]

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

    if (is.null(p_B_AD)){

   		yi <- substr(y, start, start)

      str_from <- yi[p[, 'from']] %>% 
        factor(alphabets) %>% 
        as.numeric()

      str_to <- yi[p[, 'to']] %>% 
        factor(alphabets) %>% 
        as.numeric()

      log_prob <- log_prob + log(p_AB_D[cbind(str_from, str_to, p$distance)])

    }else{

      if (start == 1){

     		yi <- substr(y, start, start)

        str_from <- yi[p[, 'from']] %>% 
          factor(alphabets) %>% 
          as.numeric()

        str_to <- yi[p[, 'to']] %>% 
          factor(alphabets) %>% 
          as.numeric()

        log_prob <- log_prob + log(p_AB_D[cbind(str_from, str_to, p$distance)])

      }else{

        yi <- substr(y, start - 1, start)

        g <- data.frame(
          i = p[, 'from'],
          j = p[, 'to'],
          from = yi[p[, 'from']],
          to = yi[p[, 'to']],
          stringsAsFactors = FALSE
        ) %>% 
          mutate(
            from2 = paste0(substr(from, 1, 1), substr(to, 1, 1)),
            to2 = paste0(substr(from, 2, 2), substr(to, 2, 2))
          ) %>%
          mutate(from = from2, to = to2) %>%
          select(-from2, -to2) 

       		str_from <- g[, 'from'] %>% 
            factor(dimers) %>% 
            as.numeric()

          str_to <- g[, 'to'] %>%
            factor(dimers) %>% 
            as.numeric()

          log_prob <- log_prob + log(p_B_AD[cbind(str_from, str_to, p$distance)])
      }
    }
  }

  log_prob <- log_prob + log(prior[p$distance])

  P <- matrix(log_prob, nrow = length(x) * (length(x) - 1) / 2, ncol = max_distance)
  post <- exp(P - apply(P, 1, logSumExp))

  p1 <- p %>% filter(distance == 1)
  D <- sparseMatrix(
    i = p1[, 'from'], 
    j = p1[, 'to'], 
    x = colSums(t(post) * (1:max_distance)), 
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
    replace_na(list(n = 0.1))  %>%
    group_by(distance, from) %>%
    mutate(prob = n / sum(n)) %>%
    select(from, to, distance, prob)

  array(
    d$prob, 
    dim = c(length(x@dimers), length(x@dimers), x@max_distance),
    dimnames = list(x@dimers, x@dimers, NULL)
  )

} # get_conditional_transition_probability


get_conditional_replacement_probability <- function(x){

  # conditional probability of single nucleotide replacement
  d <- x@df %>%
    ungroup() %>%
    filter(k == 1) %>%
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
    replace_na(list(n = 0.1))  %>%
    group_by(distance) %>%
    mutate(prob = n / sum(n)) %>%
    select(from, to, distance, prob)

  array(
    d$prob, 
    dim = c(length(x@alphabets), length(x@alphabets), x@max_distance),
    dimnames = list(x@alphabets, x@alphabets, NULL)
  )

} # get_conditional_replacement_probability
