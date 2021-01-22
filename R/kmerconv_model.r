#' KmerConvModel
#'
KmerConvModel <- function(
	n_samples = NULL,
	resblocks = 4L,
	filters = 32L,
	kernel_size = 3L,
	rate = 0.1,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		self$n_samples <- n_samples
		self$filters <- filters
		self$kernel_size <- kernel_size
		self$n_layers <- length(filters)

		self$dense_1 <- tf$keras$layers$Dense(self$filters)

		self$resblock_layers <- lapply(seq_len(resblocks), function(i) tf$keras$Sequential(list(
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Activation(activation = 'relu'),
			tf$keras$layers$Conv2D(filters = filters, kernel_size = kernel_size, strides = 1L, padding = 'same'),
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Activation(activation = 'relu'),
			tf$keras$layers$Conv2D(filters = filters, kernel_size = kernel_size, strides = 1L, padding = 'same')
		)))

		self$final <- tf$keras$layers$Dense(1L, activation = 'relu')

		function(x, training = TRUE, mask = NULL){ 

			x <- x %>% 
				self$dense_1()

			for (i in seq_len(length(self$resblock_layers))){
				d <- x %>%
					self$resblock_layers[[i - 1]]() 
				x <- x + d
			}
			x <- x %>% self$final()
			x
		}
	})
}


#' prepare_data
#'
#' @export
#'
setMethod(
	'prepare_data',
	signature(
		model = 'KmerConvModel',
		x = 'kmer_summary' 
	),
	function(
		model,
		x,
		n_sims = 10L,
		n_samples_per_sim = 100L
	){

		D <- get_replacement_distance(x)

		input <- list()
		output <- list()

		for (i in seq_len(n_sims)){

			sim <- simulate(x@config, n_samples = n_samples_per_sim)

			is_leaf <- degree(sim@graph, mode = 'out') == 0
			leaves <- names(sim@x)[is_leaf]

			xb <- prepare_input(sim@x[is_leaf], distance = D)

			db <- distances(sim@graph, v = leaves, to = leaves) %>%
				as.matrix() %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(0L) %>%
				tf$expand_dims(3L)

			input[[i]] <- xb
			output[[i]] <- db
		}

		list(
			input = input %>% tf$concat(axis = 0L),
			distance = output %>% tf$concat(axis = 0L)
		)
	}
)

#' fit
#'
#' @export
#'
setMethod(
	'fit',
	signature(
		model = 'KmerConvModel',
		x = 'tf_dataset' 
	),
	function(
		model,
		x,
		batch_size = 128L,
		epochs = 10000L,
		evaluation_data = NULL,
		test_size = 0.15,
		compile = FALSE
	){

    x <- x %>%
      split_dataset(batch_size = batch_size, test_size = test_size)

		optimizer <- tf$keras$optimizers$Adam(1e-4, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		mse <- tf$keras$losses$MeanSquaredError(reduction = 'none')

		train_step <- function(x, y){

			w <- tf$linalg$diag(tf$zeros(x$shape[[2]]), padding_value = 1)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				y_pred <- model@model(x, training = TRUE)
				loss <- (w * mse(y, y_pred)) %>% 
					tf$reduce_mean(shape(1L, 2L)) %>%
					tf$reduce_mean()
			})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss = loss
			)
		}

		test_step <- function(x, y){
			w <- tf$linalg$diag(tf$zeros(x$shape[[2]]), padding_value = 1)
			y_pred <- model@model(x, training = FALSE)
			loss <- (w * mse(y, y_pred)) %>% 
				tf$reduce_mean(shape(1L, 2L)) %>%
				tf$reduce_mean()
			list(
				loss = loss
			)
		}

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
			test_step <- tf_function(test_step) # convert to graph mode
		}

		for (epoch in 1:epochs){

			loss_train <- NULL
			loss_test <- NULL

			iter <- x$train %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()

			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$input, batch$distance)
				loss_train <- c(loss_train, as.numeric(res$loss))
			})

			iter <- x$test %>%
				make_iterator_one_shot()

			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$input, batch$distance)
				loss_test <- c(loss_test, as.numeric(res$loss))
			})

			sprintf('epoch=%6.d/%6.d | loss_train=%13.7f | loss_test=%13.7f', epoch, epochs,  mean(loss_train), mean(loss_test)) %>%
				message()

			if (!is.null(evaluation_data) && epoch %% 10 == 0){
				for (i in 1:length(evaluation_data)){
					d <- dist_sim(evaluation_data[[i]]$sequence, model, evaluation_data[[i]]$kmer_summary)
					sprintf('epoch=%6.d/%6.d | RF=%13.7f', epoch, epochs, d %>% NJ() %>% RF.dist(evaluation_data[[i]]$tree, normalize = TRUE)) %>% 
						message()
				}
			}
		}

		model
	}
)


setMethod(
	'dist_sim',
	signature(
		x = 'phyDat',
		model = 'KmerConvModel'
	),
	function(x, model, ks, batch_size = 512L, ...){

		v <- names(x)
		D <- get_replacement_distance(ks)
		x <- prepare_input(x, distance = D)
		y <- model@model(x, training = FALSE) %>% 
			tf$squeeze(shape(0L, 3L)) %>%
			as.matrix()
		rownames(y) <- colnames(y) <- v
		y
	}
)

prepare_input <- function(x, distance){

	x <- as.character(x)

	p <- expand.grid(
		from = 1:nrow(x),
		to = 1:nrow(x)
	) %>% 
		slice(rep(1:n(), each = ncol(x))) %>%
		cbind(position = rep(1:ncol(x), nrow(x) * nrow(x)))

	from <- x[cbind(p$from, p$position)]
	to <- x[cbind(p$to, p$position)]
	d <- distance[cbind(from, to)]

	d %>%
		tf$cast(tf$float32) %>%
		tf$reshape(shape(nrow(x), nrow(x), ncol(x))) %>%
		tf$expand_dims(0L)
}

