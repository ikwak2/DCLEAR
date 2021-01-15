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
			
	attention_weights <- tf$nn$softmax(scaled_attention_logits, axis = -1)  # (..., seq_len_q, seq_len_k)

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
		self$wv <- tf$keras$layers$Dense(self$d_model)
		self$wq <- tf$keras$layers$Dense(self$d_model)
		self$wk <- tf$keras$layers$Dense(self$d_model)
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

		self$embedding <- tf$keras$layers$Embedding(alphabets, filters)
		self$dropout_1 <- tf$keras$layers$Dropout(rate)
		self$pos_encoding <- tf$cast(positional_encoding(self$n_targets, self$filters), tf$float32)

		self$enc_layers <- lapply(seq_len(encoders), function(i) TransformerEncoderLayer(self$filters, num_heads = num_heads, dff = dff, rate = rate))

		self$flatten_1 <- tf$keras$layers$Flatten()

		self$dense_2 <- tf$keras$layers$Dense(1L)

		function(x1, x2, ..., training = TRUE){

			x1 <- x1 %>% self$embedding()
			x2 <- x2 %>% self$embedding()

			x <- (x1 - x2)^2
			x <- x * tf$math$sqrt(tf$cast(self$filters, tf$float32))
			x <- x + self$pos_encoding
			x <- x %>%
				self$dropout_1()

			for (i in seq_len(length(self$enc_layers))){
				x <- self$enc_layers[[i - 1]](x)
			}

			x <- x %>% 
				self$flatten_1() %>%
				self$dense_2()
			x
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
		x = 'phyDat' 
	),
	function(
		model,
		x,
		y = NULL,
		division = 16L,
		mutation_prob,
		batch_size = 128L,
		epochs = 10000L,
		steps_per_epoch = 10L,
		compile = FALSE
	){

		DEFAULT <- '0'
		DELETION <- '-'
		DROPOUT <- '*'

		alphabets <- levels(x)
		freq <- x %>%
			as.character() %>%
			factor(alphabets) %>%
			table()

		if (any(DROPOUT %in% alphabets))
			dropout_prob <- freq[DROPOUT] / sum(freq) # dropout probability
		else
			dropout_prob <- 0

		if (any(names(freq) %in% DELETION))
			freq[DELETION] <- 0

		if (any(names(freq) %in% DELETION))
			deletion <- TRUE

		if (any(names(freq) %in% DROPOUT))
			freq[DROPOUT] <- 0

		if (missing(mutation_prob))
			mutation_prob <- 1 - (freq[DEFAULT] / sum(freq))^(1 / division)

		freq[DEFAULT] <- 0
		outcome_prob <- as.numeric(freq / sum(freq))
		names(outcome_prob) <- names(freq)
		alphabets <- names(freq)

		n_targets <- x %>%
			as.character() %>%
			ncol()

		learning_rate <- tf$keras$optimizers$schedules$PolynomialDecay(
			initial_learning_rate = 1e-4,
			decay_steps = steps_per_epoch * epochs,
			end_learning_rate = 0
		)
		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		mse <- tf$keras$losses$MeanSquaredError()
		train_step <- function(x1, x2, y){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				y_pred <- model@model(x1, x2)
				loss <- mse(y, y_pred)
			})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss = loss
			)
		}

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}
		
		for (epoch in 1:epochs){

			sim <- simulate(
				n_samples = 100L, 
				n_targets = n_targets, 
				mutation_prob = mutation_prob, 
				division = division, 
				outcome_prob = outcome_prob, 
				alphabets = alphabets, 
				deletion = deletion, 
				dropout_prob = dropout_prob
			)

			xb <- sim@x %>%
			  as.character() %>%
			  factor(levels(sim@x)) %>%
				as.numeric() %>%
				matrix(nrow = length(sim@x), ncol = ncol(as.character(sim@x)), dimnames = list(names(sim@x), NULL))
			xb  <- xb  - 1

			for (s in seq_len(steps_per_epoch)){

				i1 <- sample(rownames(xb), batch_size, replace = TRUE)	# sampling some nodes
				i2 <- sample(rownames(xb), batch_size, replace = TRUE)	# sampling some nodes

				yb <- distances(sim@graph)[cbind(i1, i2)] %>%
					tf$cast(tf$float32) %>%
					tf$expand_dims(1L)

				x1 <- xb[i1, ] %>% 
					tf$cast(tf$int64) 

				x2 <- xb[i2, ] %>% 
					tf$cast(tf$int64) 

				res <- train_step(x1, x2, yb)

				sprintf('epoch=%6.d/%6.d | step=%3.d/%3.d | loss=%15.7f', epoch, epochs, s, steps_per_epoch, res$loss) %>%
					message()
			}

			if (epoch %% 10 == 0){
				d <- dist_sim(x, model)
				d %>% NJ() %>% RF.dist(y, normalize = TRUE) %>% message()
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

		alphabets <- levels(x)
		v <- names(x)

		x <- x %>%
		  as.character() %>%
		  factor(alphabets) %>%
			as.numeric() %>%
			matrix(nrow = length(x), ncol = ncol(as.character(x)))
		x <- x - 1

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

			d[[i]] <- model@model(x1, x2) %>% as.matrix() %>% c()
		}

		d <- unlist(d)
		d <- sparseMatrix(i = param[, 1], j = param[, 2], x = d, dims = c(length(v), length(v)), dimnames = list(v, v)) %>%
			as.matrix()
		d <- d + t(d)
		d
	}
)

