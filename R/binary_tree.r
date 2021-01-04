#' perfect_binary_tree
#'
#' Generate a perfect binary tree
#'
generate_perfect_binary_tree <- function(
	division = 5L
){
	num_nodes <- 2^(division + 1) - 1

	g <- lapply(1:division, function(d){
		leaves <- (2^(d)):(2^(d+ 1) - 1)
		parents <- leaves %/% 2
		cbind(leaves, parents)
	})

	g <- do.call('rbind', g)
	g <- sparseMatrix(i = g[, 1], j = g[, 2], dims = c(num_nodes, num_nodes))
	g
} # perfect_binary_tree

#' node_code
#'
node_code <- function(x){

	a <- x
	y <- x
	for (i in 1:get_division(x)){
		a <- x %*% a
		y <- y + a
	}
	y <- y + Diagonal(n = nrow(x))
	y
}

#' get_division
#'
get_division <- function(x){
	log2(nrow(x) + 1) - 1L
} # get_division

#' is_leaf
#'
is_leaf <- function(x){
	division <- get_division(x)
	1:nrow(x) %in% (2^(division)):(2^(division + 1) - 1)
}
