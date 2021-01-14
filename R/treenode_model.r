#' SimDistModel
#'
SimDistModel <- function(
	latent_dim = 10L,
	rate = 0.1,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		self$latent_dim <- latent_dim

		function(x1, x2, ..., training = TRUE){
			browser()
		}
	})
}

#' fit
#'
#' @export
#'
setMethod(
	'fit',
	signature(
		model = 'SimDistModel',
		x = NULL,
	),
	function(
		model,
		x,
		batch_size = 128L,
		epochs = 100L,
		learning_rate = 1e-3,
		compile = FALSE
	){

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		mse <- tf$keras$losses$MeanSquaredError(reduction = 'none')
		bce <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')

		browser()

		train_step <- function(x1, x2){

			observed_distance <- abs(x1 - x2) %>% 
				tf$reduce_sum(1L)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				res <- model@model(x1, x2)

				loss_recon_1 <- bce(x1, res$x1) %>%
					tf$reduce_mean(shape(0L)) 

				loss_recon_2 <- bce(x2, res$x2) %>%
					tf$reduce_mean(shape(0L)) 

				loss_recon <- loss_recon_1 + loss_recon_2

				loss_distance <- mse(observed_distance, res$distance) %>%
					tf$reduce_mean()

				loss <- loss_recon + loss_distance
			})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss_distance = loss_distance,
				loss_recon = loss_recon,
				loss = loss
			)
		}

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}
		
		for (epoch in 1:epochs){

			i1 <- sample(1:nrow(code), batch_size)	# sampling some nodes
			i2 <- sample(1:nrow(code), batch_size)	# sampling some nodes

			x1 <- x[i1, ] %>% 
				tf$cast(tf$float32)

			x2 <- x[i2, ] %>% 
				tf$cast(tf$float32)

			res <- train_step(x1, x2)

			sprintf('epoch=%6.d/%6.d | distance_loss=%15.7f | recon=%15.7f | loss=%15.7f', epoch, epochs, res$loss_distance, res$loss_recon, res$loss) %>%
				message()

		}
		model
	}
)
