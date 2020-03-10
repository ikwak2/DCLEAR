setGeneric('kmer_replacement', function(x, ...) standardGeneric('kmer_replacement'))

setMethod(
  'kmer_replacement',
  signature(
    x = 'kmer_summary'
  ),
  function(x, ...){

    d <- get_distance_probability(x)

    max_distance <- x@max_distance
    kmers <- x@kmers

    D <- matrix(c(d), nrow = length(kmers) * length(kmers), ncol = max_distance)

    D <- matrix(
      rowSums(D %*% Diagonal(x = 1:max_distance)),
      nrow = length(kmers),
      ncol = length(kmers),
      dimnames = list(kmers, kmers)
    )
    D
  }
)


#' get_distance_probability
#'
#' Compute the conditional distribution of distance according to k-mer replacement
#'
get_distance_probability <- function(x){

  d <- x@df %>%
    ungroup() %>%
    select(from, to, distance, n) %>%
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
    mutate(prob  = n / sum(n))  %>%
    mutate(prob = ifelse(is.na(prob) & from == to & distance == 1, 1, prob)) %>%
    mutate(prob = ifelse(is.na(prob) & from != to & distance == x@max_distance, 1, prob)) %>%
    replace_na(list(prob = 0))

  array(
    d$prob, 
    dim = c(length(x@kmers), length(x@kmers), x@max_distance),
    dimnames = list(x@kmers, x@kmers, NULL)
  )

} # get_distance_probability



