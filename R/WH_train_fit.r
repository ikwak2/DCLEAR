#' Train weights for WH, and output distance object
#' 
#' Train weights for WH using the given data, and fit the distance matrix for a input sequence.
#'
#' @param x input data in phyDat format
#'
#' @param X a list of k number of input data, X[[1]] ... X[[k]]. The ith data have sequence information as phyDat format in X[[i]][[1]], and tree information in X[[i]][[2]] as phylo format.
#'
#' @return a dist object
#'
#' @author Il-Youp Kwak (ikwak2@cau.ac.kr)
#'
#' @export

WH_train_fit <- function(x, X){
    InfoW = WH_train(X)  
    dist_weighted_hamming(x, InfoW, dropout=FALSE)
}

