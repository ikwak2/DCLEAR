setGeneric('sample_kmer', function(x, ...) standardGeneric('sample_kmer'))

#' function for sampling k-mers in the simulated data
#'
setMethod(
  'sample_kmer',
  signature(
    x = 'lineage_tree'
  ),
  function(
    x,
    k = 1L  # k-mer,
  ){

    X <- x@x %>% as.character()

    p <- expand.grid(
      from = 1:nrow(X),
      to = 1:nrow(X),
      start = 1:(x@n_targets - k + 1)
     ) %>%
      filter(from < to)

    D <- x@graph %>% distances(
      v = names(x@x),
      to = names(x@x)
    )
    d <- D[cbind(p[, 'from'], p[, 'to'])]

    str_from <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'from'], p[, 'start'] + j - 1)]))
    str_to <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'to'], p[, 'start'] + j - 1)]))

    data.frame(
      from = str_from,
      to = str_to,
      distance = d
    ) %>%
     group_by(from, to, distance) %>%
     tally()
  }

) # sample_kmer  
