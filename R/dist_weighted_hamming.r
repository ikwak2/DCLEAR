setGeneric('dist_weighted_hamming', function(x, y, ...) standardGeneric('dist_weighted_hamming'))

setMethod(
	'dist_weighted_hamming',
	signature(
		x = 'phyDat',
		y = 'missing'
	),
	function(x, y, method, ...){

    # compute pairwise weighted hamming distance 

	}
)

setMethod(
	'dist_weighted_hamming',
	signature(
		x = 'phyDat',
		y = 'phyDat'
	),
	function(x, y, method, ...){

    # compute the weighted hamming distance between two matrix x and y

	}
)

#' compute the weighted hamming distance between two matrices
#
setMethod(
	'dist_weighted_hamming',
	signature(
		x = 'phyDat',
		y = 'matrix'
	),
	function(x, y, method, ...){

		num_states <- nlevels(x)
		y <- y %>% phyDat(type = 'USER', levels = levels(x))
		dist_weighted_hamming(x, y, method, ...)

	}
)
