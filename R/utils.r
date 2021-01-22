#' split_dataset
#'
#' Split a tfdataset object into training and testing sets
#'
#' @param x a tfdataset object
#' @param test_size The ratio of the testing set (default 0.15)
#' @param batch_size Batch size (default: 64L)
#' @return a list that include a training and a testing dataset, where both of them are
#'          tfdataset object
#'
split_dataset <- function(x, test_size = 0.15, batch_size = 64L){

	n <- as.numeric(x$cardinality()) # total number of samples
  n_train <- as.integer((1 - test_size) * n)  # train samples

	train <- x %>%
		dataset_take(n_train) %>%
		dataset_batch(batch_size)

	test <- x %>%
		dataset_skip(n_train) %>%
		dataset_batch(batch_size)

	list(train = train, test = test)

} # split_dataset


sample_root <- function(config){
	prob <- config@frequency
	prob[config@dropout_character] <- 0
	prob[config@deletion_character] <- 0
	sample(config@alphabets, config@n_targets, replace = TRUE, prob = prob)
}


#'
as_matrix <- function(x){
  x <- x %>%
    as.character() %>%
		factor(levels(x)) %>%
		as.numeric() %>%
		matrix(nrow = length(x), ncol = ncol(as.character(x)), dimnames = list(names(x), NULL))
	x <- x - 1L # to zero-based
  x
}

get_replacement_distance <- function(x){

	g <- x@df %>%
		group_by(from, to) %>%
		summarize(mean = sum(n * distance) / sum(n))

	D <- sparseMatrix(
		i = as.numeric(factor(as.character(g$from), x@kmers)),
		j = as.numeric(factor(as.character(g$to), x@kmers)),
		x = g$mean,
		dims = c(length(x@kmers), length(x@kmers)),
		dimnames = list(x@kmers, x@kmers)
	) %>%
		as.matrix()

	D[D == 0 & row(D) != col(D)] <- x@max_distance
	D[D == 0 & row(D) == col(D)] <- 0
	D
}
