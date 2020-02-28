setGeneric('dist_replacement', function(x, replacement_matrix, ...) standardGeneric('dist_replacement'))

#' dist_replacement
#'
#' compute the kmer replacement distance matrix
#
setMethod(
	'dist_replacement',
	signature(
		x = 'phyDat',
		replacement_matrix = 'dgCMatrix'
	),
	function(x, replacement_matrix, ...){
		dist_replacement_core(x, replacement_matrix, ...)
	}
) # 


setMethod(
	'dist_replacement',
	signature(
		x = 'phyDat',
		replacement_matrix = 'missing'
	),
	function(x, replacement_matrix, ...){

		replacement_matrix <- compute_replacement_matrix(x = x, ...)
		dist_replacement_core(x, replacement_matrix, ...)
	}
) # 

#' dist_replacement_core
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

