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
