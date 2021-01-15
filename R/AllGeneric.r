setGeneric('dist_replacement', function(x, kmer_summary, ...) standardGeneric('dist_replacement'))

setGeneric('dist_weighted_hamming', function(x, wVec, dropout=FALSE, ...) standardGeneric('dist_weighted_hamming'))

setGeneric('substr_kmer', function(x, ...) standardGeneric('substr_kmer'))

setGeneric('summarize_kmer', function(x, ...) standardGeneric('summarize_kmer'))

setGeneric('fit', function(model, x, ...) standardGeneric('fit'))

setGeneric('dist_sim', function(x, model, ...) standardGeneric('dist_sim'))

setGeneric('simulate', function(config, ...) standardGeneric('simulate'))
