#' get_angles
#' Adopted from https://www.tensorflow.org/tutorials/text/transformer
#'
get_angles <- function(pos, i, d_model){
	angle_rates <- 1 / (10000^ ( (2 * i %/% 2) / d_model ))
  pos %*% angle_rates
}

#'
positional_encoding <- function(position, d_model){
	angle_rads <- get_angles(
		 matrix(0:(position - 1), position, 1),
		 matrix(0:(d_model - 1), 1, d_model),
		 d_model
	)
  even <- seq(1, ncol(angle_rads), by = 2)
	angle_rads[, even] <- sin(angle_rads[, even])
	odd <- seq(2, ncol(angle_rads), by = 2)
	angle_rads[, odd] <- cos(angle_rads[, odd])
	angle_rads
} # positional_encoding

#' scaled_dot_product_attention
#'
#' Calculate the attention weights.
#'
#' 	 q, k, v must have matching leading dimensions.
#'   k, v must have matching penultimate dimension, i.e.: seq_len_k = seq_len_v.
#'	 The mask has different shapes depending on its type(padding or look ahead) 
#'	 but it must be broadcastable for addition.
#' 
#' @param q: query shape == (..., seq_len_q, depth)
#' @param k: key shape == (..., seq_len_k, depth)
#' @param v: value shape == (..., seq_len_v, depth_v)
#' @param mask: Float tensor with shape broadcastable to (..., seq_len_q, seq_len_k). Defaults to None.
#'
#' @return output, attention_weights
#'
scaled_dot_product_attention <- function(q, k, v, mask = NULL){

	matmul_qk <- tf$matmul(q, k, transpose_b = TRUE)  # (..., seq_len_q, seq_len_k)
	dk <- tf$cast(tf$shape(k) %>% tail(1), tf$float32)
	scaled_attention_logits <- matmul_qk / tf$math$sqrt(dk)

	if (!is.null(mask)){
		scaled_attention_logits <- scaled_attention_logits + (mask * -1e9)  
	}
			
	attention_weights <- tf$nn$softmax(scaled_attention_logits, axis = -1L)  # (..., seq_len_q, seq_len_k)

	output <- tf$matmul(attention_weights, v)  # (..., seq_len_q, depth_v)

	list(output = output, attention_weights = attention_weights)

} # scaled_dot_product_attention


#' MultiHeadAttention
#'
MultiHeadAttention <- function(
	d_model, 
	num_heads,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$num_heads <- num_heads
		self$d_model <- d_model

		stopifnot(self$d_model %% self$num_heads == 0)

		self$depth <- as.integer(self$d_model / self$num_heads)
		self$wq <- tf$keras$layers$Dense(self$d_model)
		self$wk <- tf$keras$layers$Dense(self$d_model)
		self$wv <- tf$keras$layers$Dense(self$d_model)
		self$dense <- tf$keras$layers$Dense(self$d_model)

		self$split_heads <- function(x, batch_size){
			# Split the last dimension into (num_heads, depth).
			#	Transpose the result such that the shape is (batch_size, num_heads, seq_len, depth)
			x <- tf$reshape(x, c(batch_size, -1L, self$num_heads, self$depth))
			tf$transpose(x, perm = c(0L, 2L, 1L, 3L))
		}

		function(inputs, mask = NULL, ...){
			batch_size <- tf$shape(inputs$q)[1]
			q <- self$wq(inputs$q)
			k <- self$wk(inputs$k)
			v <- self$wv(inputs$v)
			q <- self$split_heads(q, batch_size)  # (batch_size, num_heads, seq_len_q, depth)
			k <- self$split_heads(k, batch_size)  # (batch_size, num_heads, seq_len_q, depth)
			v <- self$split_heads(v, batch_size)  # (batch_size, num_heads, seq_len_q, depth)
			y <- scaled_dot_product_attention(q, k, v, mask)
			scaled_attention <- y$output
			scaled_attention <- tf$transpose(scaled_attention, perm = c(0L, 2L, 1L, 3L))  # (batch_size, seq_len_q, num_heads, depth)
			concat_attention <- tf$reshape(scaled_attention, c(batch_size, -1L, self$d_model))  # (batch_size, seq_len_q, d_model)
			output <- self$dense(concat_attention)  # (batch_size, seq_len_q, d_model)
			list(output = output, attention_weights = y$attention_weights)
		}
	})
}

#' point_wise_feed_forward_network
#'
point_wise_feed_forward_network <- function(d_model, dff){
	model <- tf$keras$Sequential()
	model$add(tf$keras$layers$Dense(dff, activation = 'relu')) # (batch_size, seq_len, dff)
	model$add(tf$keras$layers$Dense(d_model)) # (batch_size, seq_len, d_model)
	model
}


#' TransformerEncoderLayer
#'
TransformerEncoderLayer <- function(
	d_model, 
	num_heads, 
	dff, 
	rate = 0.1,
	name = NULL
){
	keras_model_custom(name = name, function(self){
		self$mha <- MultiHeadAttention(d_model, num_heads)
		self$ffn <- point_wise_feed_forward_network(d_model, dff)
		self$layernorm1 <- tf$keras$layers$LayerNormalization(epsilon = 1e-3)
		self$layernorm2 <- tf$keras$layers$LayerNormalization(epsilon = 1e-3)
		self$dropout1 <- tf$keras$layers$Dropout(rate)
		self$dropout2 <- tf$keras$layers$Dropout(rate)

		function(x, training = TRUE, mask = NULL){
			res <- self$mha(list(q = x, k = x, v = x), mask)  # (batch_size, input_seq_len, d_model)
			attn_output <- res$output
			attn_output <- self$dropout1(attn_output, training = training)
			out1 <- self$layernorm1(x + attn_output)  # (batch_size, input_seq_len, d_model)
			ffn_output <- self$ffn(out1)  # (batch_size, input_seq_len, d_model)
			ffn_output <- self$dropout2(ffn_output, training = training)
			out2 <- self$layernorm2(out1 + ffn_output)  # (batch_size, input_seq_len, d_model)
			out2
		}
	})
}

#' SimDistModel
#'
SimDistModel <- function(
	n_targets = NULL,
	encoders = 3L,
	alphabets = NULL,
	num_heads = 4L,
	dff = 256L,
	filters = 32L,
	rate = 0.1,
	name = NULL
){
	keras_model_custom(name = name, function(self){

		self$filters <- filters
		self$n_targets <- n_targets
		self$alphabets <- alphabets

		self$distance_token <- alphabets
		self$input_dim <- self$distance_token + 1L

		self$sequence_length <- self$n_targets + 1L

		self$embedding <- tf$keras$layers$Embedding(self$input_dim, as.integer(self$filters / 2L))
		self$dropout_1 <- tf$keras$layers$Dropout(rate)
		self$pos_encoding <- tf$cast(positional_encoding(self$n_targets + 1L, self$filters), tf$float32)
		self$bn <- tf$keras$layers$BatchNormalization()

		self$enc_layers <- lapply(seq_len(encoders), function(i) TransformerEncoderLayer(self$filters, num_heads = num_heads, dff = dff, rate = rate))

		self$dense_1 <- tf$keras$layers$Dense(self$alphabets, activation = 'softmax')
		self$dense_2 <- tf$keras$layers$Dense(1L)


		function(x, training = TRUE, mask = NULL){ 

			batch_size <- x[[1]]$shape[[1]]

			d <- self$distance_token %>% 
				tf$constant(dtype = tf$int64, shape = shape(1L, 1L)) %>%
				tf$tile(shape(batch_size, 1L))

			sequence_1 <- x[[1]]
			sequence_2 <- x[[2]]

			sequence_1 <- list(d, sequence_1) %>% tf$concat(1L)
			sequence_2 <- list(d, sequence_2) %>% tf$concat(1L)

			sequence_1 <- sequence_1 %>% self$embedding()
			sequence_2 <- sequence_2 %>% self$embedding()

			sequence_1 <- sequence_1 * tf$math$sqrt(tf$cast(self$filters, tf$float32))
			sequence_2 <- sequence_2 * tf$math$sqrt(tf$cast(self$filters, tf$float32))

			x <- list(sequence_1, sequence_2) %>% tf$concat(2L)
			x <- x + self$pos_encoding

			x <- x %>%
				self$dropout_1(training = training)

			if (!is.null(mask)){

				mask <- list(tf$zeros(shape(batch_size, 1L)), mask) %>% 
					tf$concat(1L) %>%
					tf$reshape(shape(batch_size, 1L, 1L, self$sequence_length))
			}

			for (i in seq_len(length(self$enc_layers))){
				x <- self$enc_layers[[i - 1]](x, training = training, mask = mask)
			}

			ancestor <- x[, 2:self$sequence_length, ] %>% 
				self$dense_1()

			distance <- x[, 1, ] %>% self$dense_2()

			list(
				ancestor = ancestor,
				distance = distance
			)
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
		model = 'SimDistModel',
		x = 'lineage_tree_config' 
	),
	function(
		model,
		x,
		n_sims = 10L,
		divisions = 6L,
		n_samples_per_sim = 100L
	){

		x1 <- NULL
		x2 <- NULL
		d <- NULL
		a <- NULL

		x@division <- divisions

		for (i in seq_len(n_sims)){

			root <- sample_root(x)
			x@root <- root
			sim <- simulate(x)

			xb <- as_matrix(sim@x)
			db <- distances(sim@graph)

			df <- data.frame(
				from = rep(rownames(db), ncol(db)),
				to = rep(rownames(db), each = nrow(db)),
				distance = c(db)
			) %>%
				filter(from != to)  %>%
				filter(is_leaf[from] & is_leaf[to])

			is_leaf <- degree(sim@graph, mode = 'out') == 0
			leaves <- names(sim@x)[is_leaf]

			i1 <- sample(leaves, n_samples_per_sim, replace = TRUE)	# sampling some nodes
			i2 <- sample(leaves, n_samples_per_sim, replace = TRUE)	# sampling some nodes
			
			v <- unique(c(i1, i2))
			p <- ego(sim@graph, order = length(V(sim@graph)), nodes = v, mode = 'in')
			p <- lapply(p, names)
			
			m <- sparseMatrix(
				i = as.numeric(factor(rep(v, sapply(p, length)), V(sim@graph)$name)),
				j = as.numeric(factor(unlist(p), V(sim@graph)$name)),
				dims = c(length(V(sim@graph)), length(V(sim@graph))),
				dimnames = list( V(sim@graph)$name,  V(sim@graph)$name)
			)

			ancestor <- sapply(apply(m[i1, ]& m[i2, ], 1, which), function(xx) names(xx)[length(xx)])

			d <- c(d, db[cbind(i1, i2)])
			x1 <- rbind(x1, xb[i1, ])
			x2 <- rbind(x2, xb[i2, ])
			a <- rbind(a, xb[ancestor, ])
		}

		list(
			distance = d %>% tf$cast(tf$float32), 
			sequence_1 = x1 %>% tf$cast(tf$int64),
			sequence_2 = x2 %>% tf$cast(tf$int64),
			ancestor = a %>% tf$cast(tf$int64)
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
		model = 'SimDistModel',
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

		total_steps <- as.integer(as.numeric(x$cardinality()) * (1 - test_size) * epochs / batch_size)

    x <- x %>%
      split_dataset(batch_size = batch_size, test_size = test_size)

		learning_rate <- tf$keras$optimizers$schedules$PolynomialDecay(
			initial_learning_rate = 1e-4,
			decay_steps = total_steps,
			end_learning_rate = 0
		)
		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		mse <- tf$keras$losses$MeanSquaredError()
		cce <- tf$keras$losses$CategoricalCrossentropy(reduction = tf$keras$losses$Reduction$SUM)

		train_step <- function(sequence_1, sequence_2, ancestor, distance){

			batch_size <- sequence_1$shape[[1]]

#			mask <- (tf$random$uniform(shape(batch_size, model@model$n_targets)) < 0.15) %>% 
#				tf$cast(tf$float32) 

			ancestor <-  ancestor %>%
				tf$one_hot(model@model$alphabets) 

			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(list(sequence_1, sequence_2), training = TRUE)
				loss_distance <- mse(distance, res$distance)
				loss_ancestor <- cce(ancestor , res$ancestor)
				loss <- loss_distance + loss_ancestor
			})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss_distance = loss_distance,
				loss_ancestor = loss_ancestor,
				loss = loss
			)
		}

		test_step <- function(sequence_1, sequence_2, ancestor, distance){
			ancestor <-  ancestor %>%
				tf$one_hot(model@model$alphabets) 
			res <- model@model(list(sequence_1, sequence_2), training = FALSE)
			loss_distance <- mse(distance, res$distance)
			loss_ancestor <- cce(ancestor , res$ancestor)
			loss <- loss_distance + loss_ancestor
			list(
				loss_distance = loss_distance,
				loss_ancestor = loss_ancestor,
				loss = loss
			)
		}

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
			test_step <- tf_function(test_step) # convert to graph mode
		}

		for (epoch in 1:epochs){

			loss_train <- NULL
			loss_distance_train <- NULL
			loss_ancestor_train <- NULL
			loss_test <- NULL
			loss_distance_test <- NULL
			loss_ancestor_test <- NULL

			iter <- x$train %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()

			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$sequence_1, batch$sequence_2, batch$ancestor, batch$distance)
				loss_train <- c(loss_train, as.numeric(res$loss))
				loss_distance_train <- c(loss_distance_train, as.numeric(res$loss_distance))
				loss_ancestor_train <- c(loss_ancestor_train, as.numeric(res$loss_ancestor))
			})

			iter <- x$test %>%
				make_iterator_one_shot()

			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$sequence_1, batch$sequence_2, batch$ancestor, batch$distance)
				loss_test <- c(loss_test, as.numeric(res$loss))
				loss_distance_test <- c(loss_distance_test, as.numeric(res$loss_distance))
				loss_ancestor_test <- c(loss_ancestor_test, as.numeric(res$loss_ancestor))
			})

			sprintf('epoch=%6.d/%6.d | loss_distance_train=%13.7f | loss_ancestor_train=%13.7f | loss_train=%13.7f | loss_distance_test=%13.7f | loss_ancestor_test=%13.7f | loss_test=%13.7f', 
				epoch, epochs, 
				mean(loss_distance_train), mean(loss_ancestor_train), mean(loss_train), 
				mean(loss_distance_test), mean(loss_ancestor_test), mean(loss_test)
			) %>%
				message()

			if (!is.null(evaluation_data) && epoch %% 10 == 0){
				for (i in 1:length(evaluation_data)){
					d <- dist_sim(evaluation_data[[i]]$sequence, model)
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
		model = 'SimDistModel'
	),
	function(x, model, batch_size = 512L, ...){

		v <- names(x)

		x <- x %>% 
			as_matrix()

		param <- expand.grid(i = seq_len(nrow(x)), j = seq_len(nrow(x))) %>%
			as.matrix()
		param <- param[param[, 1] < param[, 2], ]

		starts <- seq(1, nrow(param), by = batch_size)
		ends <- starts + batch_size - 1
		ends[length(ends)] <- min(nrow(param), ends[length(ends)])

		d <- list()
		for (i in 1:length(starts)){

			b <- starts[i]:ends[i]

			x1 <- x[param[b, 1], , drop = FALSE]%>% 
				tf$cast(tf$int64) 

			x2 <- x[param[b, 2], , drop = FALSE]%>% 
				tf$cast(tf$int64) 

			res <- model@model(list(x1, x2), training = FALSE)

			d[[i]] <- res$distance %>% as.matrix() %>% c()
		}

		d <- unlist(d)
		d <- sparseMatrix(i = param[, 1], j = param[, 2], x = d, dims = c(length(v), length(v)), dimnames = list(v, v)) %>%
			as.matrix()
		d <- d + t(d)
		d
	}
)

#'
as_matrix <- function(x){
	x <- x %>%
		as.character() %>%
		factor(levels(x)) %>%
		as.numeric() %>%
		matrix(nrow = length(x), ncol = ncol(as.character(x)), dimnames = list(names(x), NULL))
	x <- x - 1L	# to zero-based
	x
}

