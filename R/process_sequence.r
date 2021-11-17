#' Process sequences
#' 
#' @param x input data in phyDat format
#' @param division cell division
#' @param dropout_character Dropout character (default: '*')
#' @param default_character Default character (default: '0')
#' @param deletion_character Deletion character (default: '-')
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
		division = 16L,
		dropout_character = '*',
		default_character = '0',
		deletion_character = '-'
	){

		lv <-  c(dropout_character, default_character, deletion_character, levels(x)) %>% unique()
		levels(x) <- lv

		alphabets <- levels(x)

		freq <- x %>%
			as.character() %>%
			factor(alphabets) %>%
			table() 

		used <- freq > 0 | names(freq) %in% c(dropout_character, default_character, deletion_character)
		alphabets <- alphabets[used]

		freq <- freq[used]

		freq <- as.numeric(freq)
		names(freq) <- alphabets

		dropout_prob <- freq[dropout_character] / sum(freq) # dropout probability

		if (freq[deletion_character] > 0){
			deletion <- TRUE
		}else
			deletion <- FALSE


		mutation_prob <- 1 - (freq[default_character] / sum(freq))^(1 / division)
		outcome_prob <- freq
		outcome_prob[deletion_character] <- 0
		outcome_prob[dropout_character] <- 0
		outcome_prob[default_character] <- 0
		outcome_prob <- outcome_prob / sum(outcome_prob)

		n_targets <- x %>%
			as.character() %>%
			ncol()

		config <- new('lineage_tree_config', 
			alphabets = alphabets,
			frequency = freq,
			mutation_prob = mutation_prob,
			outcome_prob = outcome_prob,
			n_targets = n_targets,
			dropout_prob = dropout_prob,
			division = division,
			deletion = deletion
		)

		config
	}
) # process_sequence

