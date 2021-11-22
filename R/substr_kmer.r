#' Subseting a kmer_summary object
#'
#' Summarize the short k-mer summary from the long k-mer summary
#' @param x a kmer_summary object
#' @param k k-mer length(default: 2)
#' @return a new kmer_summary object
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
  'substr_kmer',
  signature(
    x = 'kmer_summary'
  ),
  function(
    x,
    k = 2
  ){

    if (k > x@k)
      stop(sprintf('k(%d) must be smaller than the input x@k(%d)', k, x@k))
    else if (k == x@k){
      return(x)
    }else{

			from <- do.call('rbind', strsplit(x@df$from, ' '))
			to <- do.call('rbind', strsplit(x@df$to, ' '))

			x@df <- do.call('rbind', lapply(seq_len(x@k - k + 1), function(start){
				data.frame(
					from = do.call('paste', as.data.frame(from[, start:(start + k - 1), drop = FALSE])),
					to = do.call('paste', as.data.frame(to[, start:(start + k - 1), drop = FALSE])),
					distance = x@df$distance,
					n = x@df$n
				)
			})) %>%
				group_by(.data$from, .data$to, .data$distance) %>%
				summarize(n = sum(n))

			x@kmers <-   do.call('paste', do.call('expand.grid', lapply(1:k, function(j) x@alphabets)))
			x@k <- k
			x

    }
  }
)

