#' SequenceGenerator
#'
SequenceGenerator <- function(
	n_targets = NULL,
	n_alphabet = NULL,
	mlp = 128L,
	resblocks = 3L,
	filters = 32L,
	kernel_size = 3L,
	rate = 0.1,
	r = 0.3,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		self$n_targets <- as.integer(n_targets)
		self$n_alphabet <- as.integer(n_alphabet)
		self$filters <- filters

		self$mlp_layers <- lapply(seq_len(length(mlp)), function(i) tf$keras$Sequential(list(
			tf$keras$layers$Dense(units = mlp[i], activation = 'relu'),
			tf$keras$layers$Dropout(rate)
		)))

		self$dense_1 <- tf$keras$layers$Dense(self$n_targets * self$filters)

		self$resblock_layers <- lapply(seq_len(resblocks), function(i) tf$keras$Sequential(list(
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Activation(activation = 'relu'),
			tf$keras$layers$Conv1D(filters = filters, kernel_size = kernel_size, strides = 1L, padding = 'same'),
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Activation(activation = 'relu'),
			tf$keras$layers$Conv1D(filters = filters, kernel_size = kernel_size, strides = 1L, padding = 'same')
		)))

		self$dense_final <- tf$keras$layers$Dense(self$n_alphabet, activation = 'softmax')

		function(x, ..., training = TRUE){

			for (i in seq_len(length(self$mlp_layers))){
				x <- x %>%
					self$mlp_layers[[i - 1]]()
			}

			x <- x %>% 
				self$dense_1() %>%
				tf$reshape(shape(-1L, self$n_targets, self$filters))

			for (i in seq_len(length(self$resblock_layers))){
				d <- x %>%
					self$resblock_layers[[i - 1]]() 
				x <- x + r * d
			}

			x <- x %>%
				self$dense_final()
			x 	
		}
	})
}

#' TreeNodeGenerator
#'
TreeNodeGenerator <- function(
	n_nodes = NULL,
	filters = c(32L, 32L, 32L),
	kernel_size = c(3L, 3L, 3L),
	strides = c(2L, 2L, 2L),
	mlp = 32L,
	rate = 0.1,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		self$conv_layers <- lapply(seq_len(length(filters)), function(i) tf$keras$Sequential(list(
			tf$keras$layers$Conv1D(filters = filters[i], kernel_size = kernel_size[i], strides = strides[i], padding = 'valid',  activation = 'relu'),
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Dropout(rate)
		)))

		self$flatten_1 <- tf$keras$layers$Flatten()

		self$mlp_layers <- lapply(seq_len(length(mlp)), function(i) tf$keras$Sequential(list(
			tf$keras$layers$Dense(mlp[i], activation = 'relu'),
			tf$keras$layers$Dropout(rate)
		)))

		self$dense_final <- tf$keras$layers$Dense(n_nodes, activation = 'softmax')

		function(x, ..., training = TRUE){

			for (i in seq_len(length(self$conv_layers))){
				x <- x %>%
					self$conv_layers[[i - 1]]()
			}

			x <- x %>%
				self$flatten_1()

			for (i in seq_len(length(self$mlp_layers))){
				x <- x %>%
					self$mlp_layers[[i - 1]]()
			}
			x <- x %>% 
				self$dense_final()
			x	
		}

	})
}


#' SequenceDiscriminatorModel
#'
SequenceDiscriminatorModel <- function(
	resblocks = 3L,
	filters = 32L,
	kernel_size = 10L,
	rate = 0.1,
	r = 1,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		self$conv_1 <- tf$keras$Sequential(list(
			tf$keras$layers$BatchNormalization(),
			 tf$keras$layers$Activation(activation = 'relu'),
			 tf$keras$layers$Conv1D(filters = filters, kernel_size = kernel_size, strides = 1L, padding = 'same')
		))

		self$resblock_layers <- lapply(seq_len(resblocks), function(i) tf$keras$Sequential(list(
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Activation(activation = 'relu'),
			tf$keras$layers$Conv1D(filters = filters, kernel_size = kernel_size, strides = 1L, padding = 'same'),
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Activation(activation = 'relu'),
			tf$keras$layers$Conv1D(filters = filters, kernel_size = kernel_size, strides = 1L, padding = 'same')
		)))

		self$flatten_1 <-  tf$keras$layers$Flatten()
		self$dense_final <-  tf$keras$layers$Dense(1L)

		function(x, ..., training = TRUE){

			x <- x %>%
				self$conv_1()

			for (i in seq_len(length(self$resblock_layers))){
				d <- x %>%
					self$resblock_layers[[i - 1]]()
				x <- x + r * d
			}

			x %>%
				self$flatten_1() %>%
				self$dense_final()
		}
	})
}

#' TreeNodeDiscriminatorModel
#'
TreeNodeDiscriminatorModel <- function(
	mlp = c(256L, 128L, 64L),
	rate = 0.1,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		self$mlp_layers <- lapply(seq_len(length(mlp)), function(i) tf$keras$Sequential(list(
			tf$keras$layers$Dense(units = mlp[i], activation = 'relu'),
			tf$keras$layers$Dropout(rate)
		)))	
		self$dense_final <- tf$keras$layers$Dense(units = 1L)

		function(x, ..., training = TRUE){
			for (i in seq_len(length(self$mlp_layers))){
				x <- x %>%
					self$mlp_layers[[i - 1]]() 
			}
			x <- x %>%
				self$dense_final()
		}
	})
}


#' Seq2TreeModel
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
Seq2TreeModel <- function(
	latent_dim = 10L,
	n_targets = NULL,
	n_alphabet = NULL,
	n_nodes = NULL,
	rate = 0.1,
	conv_layers = 3L,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		self$n_targets <- as.integer(n_targets)
		self$n_alphabet <- as.integer(n_alphabet)

		self$sequence_generator <- SequenceGenerator(
			n_targets = n_targets,
			n_alphabet = n_alphabet,
			resblocks = 2L,
			filters = 32L,
			kernel_size = 2L,
			rate = 0.3,
			r = 1
		)

		self$node_generator <- TreeNodeGenerator(
			n_nodes = n_nodes,
			filters = c(32L, 32L, 32L),
			kernel_size = c(3L, 3L, 3L),
			strides = c(2L, 2L, 2L),
			rate = rate
		)

		self$sequence_discriminator <- SequenceDiscriminatorModel(
			resblocks = 5L,
			filters = 32L,
			kernel_size = 2L,
			rate = 0.3,
			r = 1
		)
		self$node_discriminator <- TreeNodeDiscriminatorModel(
			mlp = c(256L, 128L, 64L),
			rate = rate
		)
	})
}

#' fit
#'
#' @export
#'
setMethod(
	'fit',
	signature(
		model = 'Seq2TreeModel',
		x = 'matrix'
	),
	function(
		model,
		x,
		code,
		batch_size = 32L,
		epochs = 100L,
		learning_rate = 1e-3,
		lambda = 10,
		compile = FALSE
	){

		optimizer_sequence_generator <- tf$keras$optimizers$Adam(learning_rate, beta_1=0.5)
		optimizer_sequence_discriminator <- tf$keras$optimizers$Adam(learning_rate, beta_1=0.5)
		optimizer_node_generator <- tf$keras$optimizers$Adam(learning_rate, beta_1=0.5)
		optimizer_node_discriminator <- tf$keras$optimizers$Adam(learning_rate, beta_1=0.5)

		cce <- tf$keras$losses$CategoricalCrossentropy()

		bce <- tf$keras$losses$BinaryCrossentropy()	# for reconstruction loss

		discriminator_loss <- function(real_output, fake_output){
			bce <- tf$keras$losses$BinaryCrossentropy(from_logits = TRUE)	
			real_loss <- bce(tf$ones_like(real_output), real_output) 
			fake_loss <- bce(tf$zeros_like(fake_output), fake_output)
			real_loss + fake_loss
		}

		generator_loss <- function(fake_output){
			bce <- tf$keras$losses$BinaryCrossentropy(from_logits = TRUE)	
			bce(tf$ones_like(fake_output), fake_output)
		}

		train_step <- function(x_real, y_real){

			with(tf$GradientTape(persistent = TRUE) %as% tape, { 

				x_fake <- y_real %>%
					model@model$sequence_generator()

				y_fake <- x_real %>%
					model@model$node_generator()

				x_cycled <- y_fake %>%
					model@model$sequence_generator()

				y_cycled <- x_fake %>%
					model@model$node_generator()

				output_x_real <- model@model$sequence_discriminator(x_real)
				output_x_fake <- model@model$sequence_discriminator(x_fake)
				output_y_real <- model@model$node_discriminator(y_real)
				output_y_fake <- model@model$node_discriminator(y_fake)

				loss_gen_x <- generator_loss(output_x_fake)
				loss_gen_y <- generator_loss(output_y_fake)

				loss_cycle_x <- cce(x_real, x_cycled)
				loss_cycle_y <- bce(y_real, y_cycled)
				loss_cycle <- loss_cycle_x + loss_cycle_y

				loss_gen_sequence <- loss_gen_x + lambda * loss_cycle
				loss_gen_node <- loss_gen_y + lambda * loss_cycle

				loss_disc_sequence <- discriminator_loss(output_x_real, output_x_fake)
				loss_disc_node <- discriminator_loss(output_y_real, output_y_fake)
			})

			v <- model@model$sequence_generator$trainable_variables
			gradients <- tape$gradient(loss_gen_sequence, v)
			list(gradients, v) %>%
				purrr::transpose() %>%
				optimizer_sequence_generator$apply_gradients()

			v <- model@model$node_generator$trainable_variables
			gradients <- tape$gradient(loss_gen_node, v)
			list(gradients, v) %>%
				purrr::transpose() %>%
				optimizer_node_generator$apply_gradients()

			v <- model@model$sequence_discriminator$trainable_variables
			gradients <- tape$gradient(loss_disc_sequence, v)
			list(gradients, v) %>%
				purrr::transpose() %>%
				optimizer_sequence_discriminator$apply_gradients()

			v <- model@model$node_discriminator$trainable_variables
			gradients <- tape$gradient(loss_disc_node, v)
			list(gradients, v) %>%
				purrr::transpose() %>%
				optimizer_node_discriminator$apply_gradients()

			list(
				x_fake = x_fake,
				loss_gen_sequence = loss_gen_sequence,
				loss_gen_node = loss_gen_node,
				loss_disc_sequence = loss_disc_sequence,
				loss_disc_node = loss_disc_node,
				loss_cycle = loss_cycle
			)
		}

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}

		x <- x %>% 
			tf$cast(tf$int64) %>%
			tf$one_hot(model@model$n_alphabet) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		for (epoch in 1:epochs){

      loss_gen_sequence <- NULL
      loss_gen_node <- NULL
      loss_disc_sequence <- NULL
      loss_disc_node <- NULL
			loss_cycle <- NULL

			iter <- x %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()

			until_out_of_range({

				xb <- iterator_get_next(iter)

				i <- sample(1:nrow(code), xb$shape[[1]], replace = TRUE)	# sampling nodes
				yb <- code[i, ] %>% 
					tf$cast(tf$float32)

				res <- train_step(xb, yb)

				loss_gen_sequence <- c(loss_gen_sequence, as.numeric(res$loss_gen_sequence))
				loss_gen_node <- c(loss_gen_node, as.numeric(res$loss_gen_node))
				loss_disc_sequence <- c(loss_disc_sequence, as.numeric(res$loss_disc_sequence))
				loss_disc_node <- c(loss_disc_node, as.numeric(res$loss_disc_node))
				loss_cycle <- c(loss_cycle, as.numeric(res$loss_cycle))
			})

			res$x_fake %>%
			  tf$math$argmax(2L) %>%
			  as.matrix() %>%
				table() %>%
				print()

#			xb %>%
#			  tf$math$argmax(2L) %>%
#			  as.matrix() %>%
#				table() %>%
#				print()

			sprintf('epoch=%6.d/%6.d | loss_gen_sequence=%15.7f | loss_gen_nodee=%15.7f| loss_disc_sequence=%15.7f | loss_disc_node=%15.7f | loss_cycle=%15.7f', epoch, epochs, mean(loss_gen_sequence), mean(loss_gen_node), mean(loss_disc_sequence), mean(loss_disc_node), mean(loss_cycle)) %>%
				message()

		}
		model
	}
)
