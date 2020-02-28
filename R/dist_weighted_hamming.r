setGeneric('dist_weighted_hamming', function(x, y, ...) standardGeneric('dist_weighted_hamming'))

setMethod(
	'dist_weighted_hamming',
	signature(
		x = 'phyDat',
		y = 'missing'
	),
	function(x, y, ...){

		num_states <- nlevels(x)
		X <- as_matrix(x)
		X <- X + 1

		ws <- x %>% 
			as.character() %>% 
			factor(levels(x)) %>% 
			table()
		ws <- ws / sum(ws)

		InfoW <- -log(ws)
		DEFAULT <- '0'
		DROPOUT <- '-'

		InfoW[DEFAULT] <- 1.5
		InfoW[DROPOUT] <- 0.8
		InfoW[15:nlevels(x)] <- 7
		InfoW[is.infinite(InfoW)] <- 0

		netD <- matrix(0, nrow(X), nrow(X))

		for(i in seq_len(nlevels(x))){

			ws0 <- InfoW
			ws0[i] <- 0

			tmpx2n <- matrix(ws0[X], nrow(X), ncol(X))

			D <- InfoW[i] * tmpx2n %*% t(X == i)
			netD <- netD + (D + t(D))
		}

		netD <- netD / max(netD)

		return(as.dist(netD) )
	}
)

setMethod(
	'dist_weighted_hamming',
	signature(
		x = 'phyDat',
		y = 'phyDat'
	),
	function(x, y, method, ...){

		num_states <- nlevels(x)
		X <- as_matrix(x)
		X <- X + 1

		freq <- x %>%
			as.character() %>%
			factor(levels(x)) %>%
			table()

		freq <- freq[freq > 0]

		if (method == 'ic'){
			w <- -log(freq / sum(freq))
		}else if (method == 'none'){
			w <- rep(1, length(freq))
			names(w) <- names(freq)
		}

		d <- matrix(0, length(x), length(y), dimnames = list(names(x), names(y)))
		for (i in seq_len(length(w))){

			xi <- x %>% as.character() == names(w)[i]
			yi <- y %>% as.character() == names(w)[i]

			d <- d + w[i] * (xi %*% t(1 - yi))
		}
		d <- d / max(d)
		d
	}
)

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
