#' Tree2SeqModel
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
Tree2SeqModel <- function(
	rate = 0.1,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		function(x, ..., training = TRUE){
		}
	})
}

#' prepare_data
#'
#' @export
#'
setMethod(
	'fit',
	signature(
		model = 'Tree2SeqModel',
		x = 'phyDat'
	),
	function(
		model,
		x,
		batch_size = 6L,
		set_size = 32L,
		division = 16L,
		epochs = 100L,
		learning_rate = 1e-3,
		compile = FALSE
	){

		x <- x %>% 
			as.character() %>% 
			factor(levels(x)) %>% 
			as.numeric() %>% 
			matrix(nrow = length(x), ncol = ncol(as.character(x)))

		x <- x - 1

		g <- generate_perfect_binary_tree(division)
		g <- node_code(g)

		train_step <- function(x, y){
			browser()
			model@model(x, y)
		}
		
		for (epoch in 1:epochs){

			i <- sample(1:nrow(g), batch_size * set_size)
			j <- sample(1:nrow(x), batch_size * set_size, replace = TRUE)

			yb <- g[i, ] %>% 
				as.matrix() %>%
				tf$cast(tf$float32) %>%
				tf$reshape(shape(batch_size, set_size, -1L))

			xb <- x[j, ] %>%
				tf$cast(tf$int64) %>%
				tf$reshape(shape(batch_size, set_size, -1L))

			train_step(xb, yb)


			browser()

		}

		browser()
	}
)
