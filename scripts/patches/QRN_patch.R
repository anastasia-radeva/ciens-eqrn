#' Recurrent QRN fitting function
#'
#' @description
#' Used to fit a recurrent quantile regression neural network on a data sample.
#' Use the [QRN_fit_multiple()] wrapper instead, with `data_type="seq"`, for better stability using fitting restart.
#'
#' @param X Matrix of covariates, for training. Entries must be in sequential order.
#' @param Y Response variable vector to model the conditional quantile of, for training. Entries must be in sequential order.
#' @param q_level Probability level of the desired conditional quantiles to predict.
#' @param hidden_size Dimension of the hidden latent state variables in the recurrent network.
#' @param num_layers Number of recurrent layers.
#' @param rnn_type Type of recurrent architecture, can be one of `"lstm"` (default) or `"gru"`.
#' @param p_drop Probability parameter for dropout before each hidden layer for regularization during training.
#' @param learning_rate Initial learning rate for the optimizer during training of the neural network.
#' @param L2_pen L2 weight penalty parameter for regularization during training.
#' @param seq_len Data sequence length (i.e. number of past observations) used during training to predict each response quantile.
#' @param scale_features Whether to rescale each input covariates to zero mean and unit covariance before applying the network (recommended).
#' @param n_epochs Number of training epochs.
#' @param batch_size Batch size used during training.
#' @param X_valid Covariates in a validation set, or `NULL`. Entries must be in sequential order.
#' Used for monitoring validation loss during training, enabling learning-rate decay and early stopping.
#' @param Y_valid Response variable in a validation set, or `NULL`. Entries must be in sequential order.
#' Used for monitoring validation loss during training, enabling learning-rate decay and early stopping.
#' @param lr_decay Learning rate decay factor.
#' @param patience_decay Number of epochs of non-improving validation loss before a learning-rate decay is performed.
#' @param min_lr Minimum learning rate, under which no more decay is performed.
#' @param patience_stop Number of epochs of non-improving validation loss before early stopping is performed.
#' @param tol Tolerance for stopping training, in case of no significant training loss improvements.
#' @param fold_separation Index of fold separation or sequential discontinuity in the data.
#' @param warm_start_path Path of a saved network using [torch::torch_save()], to load back for a warm start.
#' @param patience_lag The validation loss is considered to be non-improving
#' if it is larger than on any of the previous `patience_lag` epochs.
#' @param optim_met DEPRECATED. Optimization algorithm to use during training. `"adam"` is the default.
#' @param seed Integer random seed for reproducibility in network weight initialization.
#' @param verbose Amount of information printed during training (0:nothing, 1:most important, 2:everything).
#' @param device (optional) A [torch::torch_device()]. Defaults to [default_device()].
#'
#' @return An QRN object of classes `c("QRN_seq", "QRN")`, containing the fitted network,
#' as well as all the relevant information for its usage in other functions.
#' @export
#' @importFrom coro loop
QRN_seq_fit <- function(X, Y, q_level, hidden_size=10, num_layers=1, rnn_type=c("lstm","gru"), p_drop=0,
                        learning_rate=1e-4, L2_pen=0, seq_len=10, stride=1, scale_features=TRUE, n_epochs=1e4, batch_size=256,
                        X_valid=NULL, Y_valid=NULL, lr_decay=1, patience_decay=n_epochs, min_lr=0, patience_stop=n_epochs,
                        tol=1e-4, fold_separation=NULL, warm_start_path=NULL, patience_lag=5, optim_met="adam",
                        seed=NULL, verbose=2, device=default_device()){
  
  if(!is.null(seed)){torch::torch_manual_seed(seed)}
  
  rnn_type <- match.arg(rnn_type)
  data_scaling <- process_features(X=X, intermediate_q_feature=FALSE,
                                   X_scaling=NULL, scale_features=scale_features)
  Xs <- data_scaling$X_scaled
  X_scaling <- data_scaling$X_scaling
  
  # Data Loader
  trainset <- mts_dataset(Y, Xs, seq_len, scale_Y=scale_features,
                          intermediate_quantiles=NULL,
                          fold_separation=fold_separation,
                          stride=stride, target="end",
                          device=device)
  
  Y_scaling <- trainset$get_Y_scaling()
  n_train <- length(trainset)
  trainloader <- torch::dataloader(trainset, batch_size=batch_size, shuffle=TRUE)
  
  # Validation dataset (if everything needed is given)
  do_validation <- (!is.null(Y_valid) & !is.null(X_valid))
  if(do_validation){
    data_scaling_ex <- process_features(X=X_valid, intermediate_q_feature=FALSE,
                                        X_scaling=X_scaling, scale_features=scale_features)
    validset <- mts_dataset(Y_valid, data_scaling_ex$X_scaled, seq_len, scale_Y=Y_scaling,
                            intermediate_quantiles=NULL,
                            fold_separation=NULL,
                            stride=stride, target="end",
                            device=device)
    n_valid <- length(validset)
    validloader <- torch::dataloader(validset, batch_size=batch_size, shuffle=TRUE)
  }
  
  # Instantiate network
  Dim_in <- ncol(Xs)+1
  if(!is.null(warm_start_path)){
    network <- torch::torch_load(warm_start_path, device=device)
    network$train()
  } else {
    network <- QRNN_RNN_net(type=rnn_type, nb_input_features=Dim_in, hidden_size=hidden_size,
                            num_layers=num_layers, dropout=p_drop)$to(device=device)
  }
  
  # Optimizer
  optimizer <- setup_optimizer_seq(network, learning_rate, L2_pen, optim_met=optim_met)
  
  # Train the network
  network$train()
  
  loss_log_train <- rep(as.double(NA), n_epochs)
  nb_stable <- 0
  if(do_validation){
    loss_log_valid <- rep(as.double(NA), n_epochs)
    nb_not_improving_val <- 0
    nb_not_improving_lr <- 0
  }
  curent_lr <- learning_rate
  for (e in 1:n_epochs) {
    loss_train <- 0
    coro::loop(for (b in trainloader) {
      # Forward pass
      net_out <- network(b[[1]])
      # Loss
      loss <- quantile_loss_tensor(net_out, b[[2]], q=q_level, return_agg="mean")
      # Check for bad initialization
      while(e<2 & is.nan(loss$item())){
        network <- QRNN_RNN_net(type=rnn_type, nb_input_features=Dim_in, hidden_size=hidden_size,
                                num_layers=num_layers, dropout=p_drop)$to(device=device)
        network$train()
        optimizer <- setup_optimizer_seq(network, learning_rate, L2_pen, optim_met=optim_met)
        net_out <- network(b[[1]])
        loss <- quantile_loss_tensor(net_out, b[[2]], q=q_level, return_agg="mean")
      }
      # zero the gradients accumulated in buffers (not overwritten in pyTorch)
      optimizer$zero_grad()
      # Backward pass: compute gradient of the loss with respect to model parameters
      loss$backward()
      # make one optimizer step
      optimizer$step()
      # store loss
      loss_train <- loss_train + (b[[2]]$size()[1] * loss / n_train)$item()
    })
    # Log loss
    loss_log_train[e] <- loss_train
    if(do_validation){
      network$eval()
      loss_valid <- 0
      coro::loop(for (b in validloader) {
        valid_out <- network(b[[1]])
        loss <- quantile_loss_tensor(valid_out, b[[2]], q=q_level, return_agg="mean")
        loss_valid <- loss_valid + (b[[2]]$size()[1] * loss / n_valid)$item()
      })
      loss_log_valid[e] <- loss_valid
      # Is the validation loss improving ?
      if(e>patience_lag){
        if(any(is.na(loss_log_valid[(e-patience_lag):e]))){
          if(is.na(loss_log_valid[e])){
            nb_not_improving_val <- nb_not_improving_val + 1
            nb_not_improving_lr <- nb_not_improving_lr + 1
            if(verbose>1){cat("NaN validation loss at epoch:", e, "\n")}
          }
        }else{
          if(loss_log_valid[e]>(min(loss_log_valid[(e-patience_lag):(e-1)])-tol)){
            nb_not_improving_val <- nb_not_improving_val + 1
            nb_not_improving_lr <- nb_not_improving_lr + 1
          }else{
            nb_not_improving_val <- 0
            nb_not_improving_lr <- 0
          }
        }
      }
      # Learning rate decay
      if(curent_lr>min_lr & nb_not_improving_lr >= patience_decay){
        optimizer <- decay_learning_rate(optimizer,lr_decay)
        curent_lr <- curent_lr*lr_decay
        nb_not_improving_lr <- 0
      }
      if(nb_not_improving_val >= patience_stop){
        if(verbose>0){
          cat("Early stopping at epoch:", e,", average train loss:", loss_log_train[e],
              ", validation loss:", loss_log_valid[e], ", lr=", curent_lr, "\n")
        }
        break
      }
      network$train()
    }
    
    # Tolerance stop
    if(e>1){
      if(abs(loss_log_train[e]-loss_log_train[e-1])<tol){
        nb_stable <- nb_stable + 1
      } else {
        nb_stable <- 0
      }
    }
    if(nb_stable >= patience_stop){
      if(verbose>0){
        cat("Early tolerence stopping at epoch:", e,", average train loss:", loss_log_train[e])
        if(do_validation){cat(", validation loss:", loss_log_valid[e])}
        cat(", lr=", curent_lr, "\n")
      }
      break
    }
    # Print progess
    if(e %% 100 == 0 || (e == 1 || e == n_epochs)){
      if(verbose>1){
        cat("Epoch:", e, "out of", n_epochs ,", average train loss:", loss_log_train[e])
        if(do_validation){cat(", validation loss:", loss_log_valid[e], ", lr=", curent_lr)}
        cat("\n")
      }
    }
  }
  
  network$eval()
  
  fit_qrn_ts <- list(fit_nn = network, interm_lvl = q_level, seq_len=seq_len, stride=stride,
                     train_loss = loss_log_train[1:e], X_scaling=X_scaling, Y_scaling=Y_scaling)
  
  if(do_validation){
    fit_qrn_ts$valid_loss <- loss_log_valid[1:e]
  }
  class(fit_qrn_ts) <- c("QRN_seq", "QRN")
  
  return(fit_qrn_ts)
}


#' Predict function for a QRN_seq fitted object
#'
#' @param fit_qrn_ts Fitted `"QRN_seq"` object.
#' @param X Matrix of covariates to predict the corresponding response's conditional quantiles.
#' @param Y Response variable vector corresponding to the rows of `X`.
#' @param q_level Optional, checks that `q_level == fit_qrn_ts$interm_lvl`.
#' @param crop_predictions Whether to crop out the fist `seq_len` observations (which are `NA`) from the returned matrix.
#' @param device (optional) A [torch::torch_device()]. Defaults to [default_device()].
#'
#' @return Matrix of size `nrow(X)` times `1`
#' (or `nrow(X)-seq_len` times `1` if `crop_predictions`)
#' containing the conditional quantile estimates of the corresponding response observations.
#' @export
#' @importFrom coro loop
QRN_seq_predict <- function(fit_qrn_ts, X, Y,
                            q_level=fit_qrn_ts$interm_lvl,
                            crop_predictions=FALSE,
                            stride = if (!is.null(fit_qrn_ts$stride)) fit_qrn_ts$stride else 1,
                            device=default_device()) {
  
  if (q_level != fit_qrn_ts$interm_lvl) stop("QRN q_level does not match in train and predict.")
  
  X_feats <- process_features(X, intermediate_q_feature=FALSE, X_scaling=fit_qrn_ts$X_scaling)$X_scaled
  
  testset <- mts_dataset(Y, X_feats, fit_qrn_ts$seq_len,
                         scale_Y=fit_qrn_ts$Y_scaling,
                         intermediate_quantiles=NULL,
                         fold_separation=NULL,
                         stride=stride, target="end",
                         device=device)
  
  bs <- min(256, length(testset))
  testloader <- torch::dataloader(testset, batch_size=bs, shuffle=FALSE)
  
  network <- fit_qrn_ts$fit_nn
  network$eval()
  
  preds <- c()
  coro::loop(for (b in testloader) {
    out <- network(b[[1]])
    preds <- c(preds, as.numeric(out))
  })
  preds <- matrix(preds, ncol=1)
  
  if (crop_predictions) {
    return(preds)
  } else {
    # full-length alignment: NA except at target indices
    out_full <- matrix(as.double(NA), nrow=nrow(X), ncol=1)
    t_idx <- testset$get_target_idx()
    out_full[t_idx, 1] <- preds[, 1]
    return(out_full)
  }
}

