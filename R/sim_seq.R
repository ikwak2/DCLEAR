#'
#' sim_seqdata 
#'
#' Generate singe cell barcode data set with tree shaped lineage information 
#'
#' @param sim_n Number of cell samples to simulate. 
#'
#' @param m Number of targets. 
#'
#' @param mu_d Mutation rate. (a scalar or a vector)
#'
#' @param d Number of cell divisions.
#'
#' @param n_s Number of possible outcome states
#'
#' @param outcome_prob Outcome probability vector (default is NULL)
#'
#' @param p_d Dropout probability
#'
#' @return The result is a list containing two objects, 'seqs' and 'tree'. The 'seqs' is 'phyDat' object of 'sim_n' number of simulated barcodes corresponding to each cell, and The 'tree' is a 'phylo' object, a ground truth tree structure for the simulated data. 
#'
#' @author Il-Youp Kwak
#'
#' @examples
#'
#' library(DCLEAR)
#' library(phangorn)
#' library(ape)
#' 
#' set.seed(1)
#' mu_d1 = c( 30, 20, 10, 5, 5, 1, 0.01, 0.001)
#' mu_d1 = mu_d1/sum(mu_d1)
#' simn = 10 # number of cell samples
#' m = 10  ## number of targets
#' sD = sim_seqdata(sim_n = simn, m = m, mu_d = 0.03,
#'         d = 12, n_s = length(mu_d1), outcome_prob = mu_d1, p_d = 0.005 )
#' ## RF score with hamming distance
#' D_hm = dist.hamming(sD$seqs)
#' tree_hm = NJ(D_hm)
#' RF.dist(tree_hm, sD$tree, normalize = TRUE)
#'
#' ## RF score with weighted hamming
#' InfoW = -log(mu_d1)
#' InfoW[1:2] = 1
#' InfoW[3:7] = 4.5
#' D_wh = dist_weighted_hamming(sD$seqs, InfoW, dropout=FALSE)
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
#' @export

sim_seqdata <- function(sim_n = 200, # number of samples to simulate 
                     m = 200,  # number of targets 
                     mu_d = 0.03,   # mutation rate (a scalar or a vector)
                     d = 15, # number of cell divisons
                     n_s = 23, # number of possible outcome states 
                     outcome_prob = NULL, # outcome probability vector
                     p_d = 0.003 # dropout prob
                     
                     ){
    
    if (length(mu_d) != m && length(mu_d) != 1  ) {
        stop("length of mu_d must be the number of targets or 1.")
    }
    
    if (!is.null(outcome_prob) ) {
        if (length(outcome_prob) != n_s ) {
        stop("length of outcome_prob must be the number of possible mutational outcomes.")
        }
    }

    acell = as.integer(rep(1,m))
    nmstrings = c( '0', '-', LETTERS[1:n_s] ) 
    sim_tree = list()
    sim_tree[[1]] = acell
    
    k = 2
    while(k < 2^d) {
        ## codes for point mutation
        mother_cell = sim_tree[[k%/%2]]
        mu_loc = runif(m) < mu_d
        mutation_cites = (mother_cell == 1) &  mu_loc
        n_mut = sum(mutation_cites)
        if (n_mut > 0) {
            if (is.null(outcome_prob)) {
                mother_cell[mutation_cites] = as.integer(sample(n_s, n_mut, replace = T)+2)
            } else {
                mother_cell[mutation_cites] = as.integer(sample(n_s, n_mut, replace = T, prob = outcome_prob)+2)
            }
        }

        ## codes for dropout
        dropout_cites = runif(m) < p_d
        if (sum(dropout_cites) > 2 ) {
            dropout_between = sample(which(dropout_cites), 2 )
            mother_cell[dropout_between[1]:dropout_between[2]] = as.integer(2)
        }

        sim_tree[[k]] = mother_cell
        k = k+1
    }
        
    sampled_leaves = sample(2^(d-1):(2^d -1) , sim_n)
    names(sampled_leaves) = sprintf('cell_%s', str_pad(1:sim_n, 4, pad = "0"))
    
    edges = matrix(0,2^d,2)
    indexes = sampled_leaves
    idx = sim_n
    iterk = sim_n
    edge_i = 1
    names(indexes) = 1:length(indexes)

    while( sum(indexes) > 1 ) {

        indexes = indexes %/% 2
        tb1 = table(indexes)

        info_tmp = which(tb1 == 2)
        if(sum(tb1 == 2) != 0) {
            for (i in 1:length(info_tmp))  {
                node_i = names(info_tmp)[i] 
                if (node_i != 0) { 

                    edges[edge_i:(edge_i+1),1] = as.integer(idx + 1)
                    edges[edge_i:(edge_i+1),2] = as.integer(names(indexes[indexes == node_i]))
                    edge_i = edge_i + 2

                    indexes[indexes == node_i] = 0
                    indexes[idx+1] = as.integer(node_i)
                    idx = idx + 1
                }
            }
            names(indexes) = 1:length(indexes)
        }
    }
    tree_ex = list()
    tree_ex$edge = edges[1:(edge_i-1),]
    tree_ex$tip.label = names(sampled_leaves)
    tree_ex$Nnode = length(unique(tree_ex$edge[,1]))
    attr(tree_ex,'class') = 'phylo'

    fixtree <- function(tree) {
        tmp_tree_nw = write.tree(tree, "")
        newtree2 = read.tree(text = tmp_tree_nw)
        return(newtree2)
    }

    tree_ex2 = fixtree(tree_ex)
    
    seqs <- sampled_leaves %>% map(~nmstrings[sim_tree[[.]]])
		seqs <- do.call('rbind', seqs)
		seqs <- seqs %>% phyDat(type = 'USER', levels = nmstrings)
    
    return(list(seqs = seqs, tree = tree_ex2))
    
}

        
