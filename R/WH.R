#'
#' WH 
#'
#' implementation of weighted hamming algorithm
#'
#' @param x Sequence object of 'phyDat' type.
#'
#' @param InfoW Weight vector for the calculation of weighted hamming distance
#'
#' @param dropout Different weighting strategy is taken to consider interval dropout with dropout = 'TRUE'. Default is, dropout = 'FALSE'.
#'
#' @return Calculated distance matrix of input sequences. The result is a 'dist' class object.  
#'
#' @author Il-Youp Kwak
#'
#' @examples
#'
#' set.seed(1)
#' library(phangorn)
#' mu_d1 = c( 30, 20, 10, 5, 5, 1, 0.01, 0.001)
#' mu_d1 = mu_d1/sum(mu_d1)
#' simn = 10 # number of cell samples
#' m = 10  ## number of targets
#' sD = sim_seqdata(sim_n = simn, m = m, mu_d = 0.03,
#'         d = 12, n_s = length(mu_d1), outcome_prob = mu_d1, p_d = 0.005 )
#'
#' ## RF score with hamming distance
#' D_h = dist.hamming(sD$seqs)
#' tree_h= NJ(D_h)
#' RF.dist(tree_h, sD$tree, normalize = TRUE)
#'
#' ## RF score with weighted hamming
#' InfoW = -log(mu_d1)
#' InfoW[1:2] = 1
#' InfoW[3:7] = 4.5
#' 
#' D_wh = WH(sD$seqs, InfoW)
#' tree_wh= NJ(D_wh)
#' RF.dist(tree_wh, sD$tree, normalize = TRUE)
#'
#' ## RF score with weighted hamming, cosidering dropout situation
#' nfoW = -log(mu_d1)
#' InfoW[1] = 1
#' InfoW[2] = 12
#' InfoW[3:7] = 3
#'
#' D_wh2 = WH(sD$seqs, InfoW, dropout=TRUE)
#' tree_wh2= NJ(D_wh2)
#' RF.dist(tree_wh2, sD$tree, normalize = TRUE)
#'
#' @export

WH <- function(x, InfoW, dropout = FALSE) {
    states <- attributes(x)$allLevels
    states2num = 1:length(states)
    names(states2num) = states
    x2 = do.call('rbind', x %>% map(~states2num[.]) )
    if (dropout) {
        D = dist_w_ham2(x2, InfoW)
    } else {
        D = dist_w_ham(x2, InfoW)
    }

    D = (D + t(D)) / max(D)
    rownames(D) = colnames(D) = names(x)
    return(as.dist(D) )
}

        
