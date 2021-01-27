#' Generic function for dist_replacement
#' @param x a sequence object
#' @param kmer_summary  a kmer_summary object
#' @param k k-mer length
#' @param ... additional parameters
#'
setGeneric('dist_replacement', function(x, kmer_summary, k, ...) standardGeneric('dist_replacement'))

#' Generic function for dist_weighted_hamming
#' @param x a sequence object
#' @param wVec weight vector
#' @param ... additional parameters
#'
setGeneric('dist_weighted_hamming', function(x, wVec, ...) standardGeneric('dist_weighted_hamming'))

#' Generic function for substr_kmer
#' @param x a kmer object
#' @param ... additional parameters
#'
setGeneric('substr_kmer', function(x, ...) standardGeneric('substr_kmer'))

#' Generic function for summarize_kmer
#' @param x a sequence object
#' @param ... additional parameters
#'
setGeneric('summarize_kmer', function(x, ...) standardGeneric('summarize_kmer'))

#' Generic function for simulate
#' @param config a lineage_tree_config object
#' @param ... additional parameters
#'
setGeneric('simulate', function(config, ...) standardGeneric('simulate'))

#' Generic function for process_sequence
#' @param x a sequence object
#' @param ... additional parameters
#'
setGeneric('process_sequence', function(x, ...) standardGeneric('process_sequence'))

#' Generic function for as_phylo
#' @param x a graph object
#' @param ... additional parameters
#'
setGeneric('as_phylo', function(x, ...) standardGeneric('as_phylo'))
