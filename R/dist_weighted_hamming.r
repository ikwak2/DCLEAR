#'
#' dist_weighted_hamming
#'
#' implementation of weighted hamming algorithm
#'
#' @param x Sequence object of 'phyDat' type.
#'
#' @param wVec Weight vector for the calculation of weighted hamming distance
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
#' mu_d1 = c( 30, 20, 10, 5, 5, 1, 0.01, 0.001)
#' mu_d1 = mu_d1/sum(mu_d1)
#' simn = 100 # number of cell samples
#' m = 200  ## number of targets
#' sD = sim_seqdata(sim_n = simn, m = m, mu_d = 0.03, d = 12, n_s = length(mu_d1), outcome_prob = mu_d1, p_d = 0.005 )
#' ## RF score with hamming distance
#' D_hm = dist.hamming(sD$seqs)
#' tree_hm = NJ(D_hm)
#' RF.dist(tree_hm, sD$tree, normalize = TRUE)
#'
#' ## RF score with weighted hamming
#' InfoW = -log(mu_d1)
#' InfoW[1:2] = 1
#' InfoW[3:7] = 4.5
#' D_wh = dist_weighted_hamming(sD$seqs, InfoW)
#' tree_wh = NJ(D_wh)
#' RF.dist(tree_wh, sD$tree, normalize = TRUE)
#'
#' ## RF score with weighted hamming, cosidering dropout situation
#' nfoW = -log(mu_d1)
#' InfoW[1] = 1
#' InfoW[2] = 12
#' InfoW[3:7] = 3
#' D_wh2 = dist_weighted_hamming(sD$seqs, InfoW, dropout = TRUE)
#' tree_wh2= NJ(D_wh2)
#' RF.dist(tree_wh2, sD$tree, normalize = TRUE)
#'
#' 
#' @export



setGeneric('dist_weighted_hamming', function(x, wVec, dropout, ...) standardGeneric('dist_weighted_hamming'))

setMethod(
	'dist_weighted_hamming',
	signature(
            x = 'phyDat',
#            y = 'missing',
            wVec = 'numeric',
            dropout = 'logical'
	),
	function(x, wVec, dropout=FALSE, ...){

            ## compute pairwise weighted hamming distance

            states <- c('0', '-', LETTERS)
            states2num = 1:length(states)
            names(states2num) = states
            x2 = do.call('rbind', x %>% map(~states2num[.]) )
            
            if (dropout) {
                D = dist_w_ham2(x2, wVec)
            } else {
                D = dist_w_ham(x2, wVec)
            }
            
            D = (D + t(D)) / max(D)
            rownames(D) = colnames(D) = names(x)
            return(as.dist(D) )

	}
)



#setMethod(
#	'dist_weighted_hamming',
#	signature(
#		x = 'phyDat',
#		y = 'phyDat'
#	),
#	function(x, y, method, ...){

    # compute the weighted hamming distance between two matrix x and y

#	}
#)

#' compute the weighted hamming distance between two matrices
#
#setMethod(
#	'dist_weighted_hamming',
#	signature(
#		x = 'phyDat',
#		y = 'matrix'
#	),
#	function(x, y, method, ...){
#
#		num_states <- nlevels(x)
#		y <- y %>% phyDat(type = 'USER', levels = levels(x))
#		dist_weighted_hamming(x, y, method, ...)
#
#	}
#)
