#' rbind
#'
#' Concatenate multiple phyDat objects
#'
#' @param ... a list of phyDat objects
#'
#' @return a phyDat object
#'
#' @export
#'
setMethod(
	'rbind',
	signature(
		... = 'phyDat'
	),
	function(
		...
	){

		# need to make sure that the levels of each phyDat are the same
		x <- list(...)
		lv <- levels(x[[1]])
		do.call('rbind', lapply(x, as.character)) %>% phyDat(type = 'USER', levels = lv)

	}
) # rbind
