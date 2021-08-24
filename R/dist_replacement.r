#' Compute the kmer replacement distance
#' 
#' Compute the kmer replacement distance between sequences
#'
#' @param x input data in phyDat format
#' @param kmer_summary a kmer_summary object
#' @param k k-mer length
#' @param ... other arguments passed to substr_kmer
#' @return a dist object
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'dist_replacement',
	signature(
		x = 'phyDat',
    kmer_summary = 'kmer_summary',
		k = 'integer'
	),
	function(x, kmer_summary, k = 2, ...){

		kmer_summary <- substr_kmer(kmer_summary, k = k, ...)

    dist_kmer_replacement_inference(x, kmer_summary, k)

  }
)


#' Compute the kmer replacement distance
#' 
#' Compute the kmer replacement distance between sequences
#'
#' @param x input data in phyDat format
#' @param kmer_summary a kmer_summary object
#' @param k k-mer length
#' @param ... other arguments passed to substr_kmer
#' @return a dist object
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'dist_replacement',
	signature(
		x = 'phyDat',
    kmer_summary = 'missing',
		k = 'integer'
	),
	function(x, kmer_summary, k = 2L, ...){

		kmer_summary <- x %>% summarize_kmer(k = k, ...)

    dist_kmer_replacement_inference(x, kmer_summary, k)

  }
)




#' Core function of computing kmer replacement distance
#' 
#' Compute the sequence distance matrix using inferred kmer replacement matrix 
#'
#' @param x input data in phyDat format
#' @param kmer_summary a kmer_summary object
#' @param k k-mers (default k=2)
#' @return a dist object
#' @author Wuming Gong (gongx030@umn.edu)
#' @importFrom rlang .data
#'
dist_kmer_replacement_inference <- function(x, kmer_summary, k = 2){

	num_alphabets <- length(kmer_summary@alphabets)
	num_kmers <- length(kmer_summary@kmers)
  sequence_length <- ncol(x %>% as.character())

  log_p_D <- log(get_distance_prior(kmer_summary))

  log_p_AB_D <- kmer_summary %>%
    substr_kmer(k = 1) %>%
    get_replacement_probability()
	log_p_AB_D <- log(log_p_AB_D + 1e-10)

	dim(log_p_AB_D) <- c(num_alphabets * num_alphabets, kmer_summary@max_distance)

	# an array for mapping alphabet to 2-mer
  km1 <- matrix(1:(num_alphabets * num_alphabets), num_alphabets, num_alphabets, dimnames = list(kmer_summary@alphabets, kmer_summary@alphabets))

  if (k > 1){
    log_p_B_AD <- get_transition_probability(kmer_summary)
		log_p_B_AD <- log(log_p_B_AD + 1e-10)
		dim(log_p_B_AD) <- c(length(kmer_summary@kmers) * length(kmer_summary@kmers), kmer_summary@max_distance)
  	km2 <- matrix(1:(num_kmers * num_kmers), num_kmers, num_kmers, dimnames = list(kmer_summary@kmers, kmer_summary@kmers))
	}

	p <- expand.grid(
		from = 1:length(x),
		to = 1:length(x)
	) %>% 
		filter(.data$from < .data$to)

  # input sequence 
	y <- do.call('rbind', strsplit(as.character(x), '')) %>%
		factor(kmer_summary@alphabets) %>% 
		as.integer() %>%
		matrix(length(x), sequence_length)

  d <- rep(0, nrow(p))

  for (start in seq_len(sequence_length - k + 1)){  # for each k-mer segment

		if (start == 1 || start %% 10 == 0){
			sprintf('posterior probability | position=%5.d/%5.d', start, sequence_length) %>% message()
		}

    log_prob <- matrix(log_p_D, nrow(p), kmer_summary@max_distance, byrow = TRUE)

    yi <- y[, start]
    str_from <- yi[p[, 'from']] 
    str_to <- yi[p[, 'to']] 

    log_prob <- log_prob + log_p_AB_D[km1[cbind(str_from, str_to)], ]

    if (k > 1){

      for (i in seq_len(k - 1)){

        start_i <- start + i - 1

        if (i > 1){
          t1 <- y[, start_i]
				  str_from <- km1[cbind(t1[p[, 'from']], t1[p[, 'to']])]
        }

        t2 <- y[, start_i + 1]
				str_to <- km1[cbind(t2[p[, 'from']], t2[p[, 'to']])]

        log_prob <- log_prob + log_p_B_AD[km2[cbind(str_from, str_to)], ]

      }
    }

    post <- exp(log_prob - rowLogSumExps(log_prob))

		di <- rowSums(post %*% diag(x = 1:kmer_summary@max_distance))

		d <- d + di

  }

	D <- sparseMatrix(
  	i = p$from,
		j = p$to,
		x = d,
		dims = c(length(x), length(x)),
		dimnames = list(names(x), names(x))
	) %>%
		as.matrix()

  D <- t(D) + D
  D %>% as.dist()

} # dist_kmer_replacement_inference

 
#' get_distance_prior 
#'
#' prior distribution of distance
#'
#' @param x a kmer_summary object
#' @return a probabilistic vector of the distribution of nodal distances
#'
#' @author Wuming Gong (gongx030@umn.edu)
#' @importFrom rlang .data
#'
get_distance_prior <- function(x){

  d <- x@df %>%
    ungroup() %>%
    group_by(.data$distance) %>%
    summarize(n = sum(n)) %>%
    mutate(prob = n / sum(n)) %>%
    select(.data$distance, .data$prob)

  prior <- rep(0, x@max_distance)
  prior[d$distance] <- d$prob
  prior

} # get_distance_prior


#' get_transition_probability
#'
#' Compute p(A,X|B,Y,d), the conditional probability of seeing a replacement from A to B
#' given the previous replacement B from Y at nodal distance d
#'
#' @param x a kmer_summary object
#' @return an 3D probabilistic array (kmers by kmers by distances)
#' @export
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
get_transition_probability <- function(x){

  # conditional transition probability
  d <- x@df %>%
    ungroup() %>%
    mutate(
      from2 = paste0(substr(.data$from, 1, 1), substr(.data$to, 1, 1)),
      to2 = paste0(substr(.data$from, 2, 2), substr(.data$to, 2, 2))
    ) %>%
    mutate(from = .data$from2, to = .data$to2) %>%
    select(.data$from, .data$to, .data$distance, .data$n) %>%
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
    mutate(n = ifelse(.data$from == .data$to, n + 1, n)) %>%
    group_by(.data$distance, .data$from) %>%
    mutate(prob = n / sum(n)) %>%
    select(.data$from, .data$to, .data$distance, .data$prob) %>%
		mutate(
			from = factor(.data$from, x@kmers),
			to = factor(.data$to, x@kmers)
		) %>%
		arrange(.data$distance, .data$to, .data$from)

  array(
    d$prob, 
    dim = c(length(x@kmers), length(x@kmers), x@max_distance),
    dimnames = list(x@kmers, x@kmers, NULL)
  )

} # get_transition_probability



#' get_replacement_probability
#'
#' Compute p(A,B|d), the conditional probability of seeing a replacement of from kmer A to B or vice versa
#'
#' @param x a kmer_summary object
#' @return an 3D probabilistic array (kmers by kmers by distances)
#'
#' @export
#'
#' @author Wuming Gong (gongx030@umn.edu)
#' @importFrom rlang .data
#'
get_replacement_probability <- function(x){

  d <- x@df %>% 
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
    mutate(n = ifelse(.data$from == .data$to, n + 1, n)) %>%
    group_by(.data$distance) %>%
    mutate(prob = n / sum(n)) %>%
    select(.data$from, .data$to, .data$distance, .data$prob) %>%
		mutate(
			from = factor(.data$from, x@kmers),
			to = factor(.data$to, x@kmers)
		) %>%
		arrange(.data$distance, .data$to, .data$from)

  array(
    d$prob, 
    dim = c(length(x@kmers), length(x@kmers), x@max_distance),
    dimnames = list(x@kmers, x@kmers, NULL)
  )

} # get_replacement_probability


