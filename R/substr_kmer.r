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
      stop(sprintf('k=%d must be smaller than the input k=%d', k, x@k))
    else if (k == x@k){
      return(x)
    }else{

      x@df <- do.call('rbind', lapply(seq_len(x@k - k + 1), function(start){
        x@df %>%
          ungroup() %>%
          mutate(
            from = substr(from, start, start + k - 1),
            to = substr(to, start, start + k - 1)
          )
      })) %>%
        group_by(from, to, distance) %>%
        summarize(n = sum(n))

      x@kmers <-   do.call('paste0', do.call('expand.grid', lapply(1:k, function(j) x@alphabets)))
      x@k <- k
      x
    }
  }
)

