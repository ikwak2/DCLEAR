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
    n_nodes = 10L,  # number of randomly sampled nodes (leaves or internval nodes)
    k = 1L  # k-mer,
  ){

    # sampling nodes
#    s <- x@x %>%
#      length() %>%
#      sample(n_nodes)

    s <- which(degree(x@graph, mode = 'out') == 0)

    X <- x@x[s] %>% as.character()

    p <- expand.grid(
      from = 1:n_nodes, 
      to = 1:n_nodes, 
      start = 1:(x@n_targets - k + 1)
     ) %>%
      filter(from < to)

    D <- x@graph %>% distances(
      v = names(x@x)[s], 
      to = names(x@x)[s]
    )
    d <- D[cbind(p[, 'from'], p[, 'to'])]
    # the k-mer at each position for each pairs of sample
    str_from <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'from'], p[, 'start'] + j - 1)]))
    str_to <- do.call('paste0', lapply(1:k, function(j) X[cbind(p[, 'to'], p[, 'start'] + j - 1)]))

    data.frame(
      from = str_from,
      to = str_to,
      distance = d
    )
  }

) # sample_kmer  
