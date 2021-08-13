#' Process sequences
#' 
#' @param x input data in phyDat format
#' @param division cell divisions (default: 16L)
#' @return a `lineage_tree_config` object
#' @author Wuming Gong (gongx030@umn.edu)
#' @export
#'
setMethod(
	'process_sequence',
	signature(
		x = 'phyDat'
	),
	function(
		x, 
		division = 16L
	){

		config <- new('lineage_tree_config')

		if (!all(levels(x) %in% config@alphabets)){
			stop('unknown characters in the input sequences')
		}

		freq <- x %>%
			as.character() %>%
			factor(config@alphabets) %>%
			table() 

		used <- freq > 0 | names(freq) %in% c(config@dropout_character, config@default_character, config@deletion_character)
		config@alphabets <- config@alphabets[used]
		config@outcome_prob <- config@outcome_prob[used]
		freq <- freq[used]

		freq <- as.numeric(freq)
		names(freq) <- config@alphabets

		dropout_prob <- freq[config@dropout_character] / sum(freq) # dropout probability

		if (freq[config@deletion_character] > 0){
			deletion <- TRUE
		}else
			deletion <- FALSE

		mutation_prob <- 1 - (freq[config@default_character] / sum(freq))^(1 / division)

		outcome_prob <- freq
		outcome_prob[config@deletion_character] <- 0
		outcome_prob[config@dropout_character] <- 0
		outcome_prob[config@default_character] <- 0
		outcome_prob <- outcome_prob / sum(outcome_prob)

		n_targets <- x %>%
			as.character() %>%
			ncol()

		config@frequency <- freq
		config@mutation_prob <- mutation_prob
		config@outcome_prob <- outcome_prob
		config@n_targets <- n_targets
		config@dropout_prob <- dropout_prob
		config@division <- division
		config@deletion <- deletion

		config
	}
) # process_sequence
