setGeneric('as_matrix', function(x, ...) standardGeneric('as_matrix'))

setMethod(
	'as_matrix',
	signature(
		x = 'phyDat'
	),
	function(x, ...){
		states <- levels(x)
		num_states <- nlevels(x)
		X <- as.character(x)
		X <- X %>%
			factor(states) %>%
			as.numeric() %>%
			matrix(nrow = nrow(X), ncol = ncol(X), dimnames = list(names(x), NULL))

		# conver the levels to zero-based so that the one hot encoding (k_one_hot) for the first level (e.g. '0') is [1,0,0, ...]
		# otherwise it will be [0, 1, 0, ...]
		X <- X - 1  
		X
	}
)


