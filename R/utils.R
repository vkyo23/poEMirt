#' @description Imputing function
#' @importFrom stats median
#' @keywords internal
#' @noRd
med_impute <- function(x) {
  x[is.na(x)] <- stats::median(x, na.rm = TRUE)
  return(x)
}

#' @description Array-handling function
#' @keywords internal
#' @noRd
dimstat <- function(array, fun, ...) apply(array, c(1, 2), fun, ...)

#' @description Generate initial values
#' @importFrom stats qlogis prcomp
#' @importFrom Rcpp sourceCpp
#' @useDynLib poEMirt, .registration = TRUE
#' @keywords internal
#' @noRd
make_init <- function(data, priors, constraint = NULL) {
  # alpha and beta
  alpha_init <- beta_init <- matrix(0, data$size$J, data$size$maxK-1)
  sums <- colSums(data$data$response, na.rm = TRUE)
  for (j in 1:data$size$J) {
    #Kj <- max(data$categories[[j]])
    cat <- data$categories[[j]]+1
    prob <- sums[j, ] / sum(sums[j, ])
    prob <- prob[cat]
    #prob <- prob[prob != 0]
    running <- cumsum(prob) / sum(prob)
    sb <- prob / (1 - (running - prob))
    sb <- sb[1:(length(sb)-1)]
    alpha_init[j, cat[1:(length(cat)-1)]] <- suppressWarnings(stats::qlogis(sb))
    beta_init[j, cat[1:(length(cat)-1)]] <- 0.1
  }
  alpha_init[is.infinite(alpha_init) | is.na(alpha_init)] <- 0.1
  
  # theta
  Yimp <- data$data$response[, , 1] / data$data$trial
  Yimp <- apply(Yimp, 2, med_impute)
  theta_init <- scale(stats::prcomp(Yimp)$x[, 1])
  if (is.null(constraint)) {
    if (exists("dynamic", data)) {
      constraint <- rep(which.max(theta_init), data$size$T)
    } 
  }
  if (exists("dynamic", data)) {
    # Static IRT
    static_fit <- poEMirtbase_fit(
      Y = data$data$response,
      S = data$data$modelinput$S, 
      Nks = data$data$modelinput$Nks, 
      alpha_old = alpha_init, 
      beta_old = beta_init, 
      theta_old = theta_init, 
      unique_categories = data$categories, 
      a0 = priors$a0, 
      A0 = priors$A0,  
      b0 = priors$b0, 
      B0 = priors$B0, 
      constraint = constraint[1] - 1, 
      std = TRUE, 
      maxit = 500, 
      verbose = 501, 
      tol = 1e-3,
      compute_ll = FALSE
    )
    theta_init <- matrix(static_fit$theta, data$size$I, data$size$T)
    alpha_init <- static_fit$alpha
    beta_init <- static_fit$beta
  } else {
    if (!is.null(constraint)) {
      if (theta_init[constraint] < 0) theta_init <- -theta_init
    }
  }
  
  L <- list(
    theta = theta_init,
    alpha = alpha_init,
    beta = beta_init
  )
  return(L)
}

#' @description Suppress messages
#' @keywords internal
#' @noRd
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#' @description Bootstrap function
#' @useDynLib poEMirt, .registration = TRUE
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @keywords internal
#' @noRd
poEMirt_boot <- function(fit, iter, verbose, save_item_parameters, thread, seed) 
{
  if (fit$info$model == "static") {
    theta_store <- matrix(NA, fit$info$data$size$I, iter)
  } else {
    theta_store <- array(NA, c(fit$info$data$size$I, fit$info$data$size$T, iter))
  }
  if (save_item_parameters) {
    alpha_store <- beta_store <- array(NA, c(fit$info$data$size$J, fit$info$data$size$maxK-1, iter))
  }
  if (thread == 1) {
    for (b in 1:iter) {
      # Simulated response
      dd <- fit$info$data
      dd$data$response <- predict.poEMirtFit(fit, type = "response", return = "array")
      dd$data$modelinput <- construct_sb_auxs(dd$data$response, dd$data$trial, dd$categories)
      fit_boot <- quiet(
        poEMirt(
          data = dd,
          model = fit$info$model,
          constraint = fit$info$constraint,
          alpha_fix = fit$info$alpha_fix,
          theta_std = fit$info$theta_std,
          priors = fit$info$priors,
          control = fit$info$control
        )
      )
      if (fit$info$model == "static") {
        theta_store[, b] <- fit_boot$parameter$theta
      } else {
        theta_store[, , b] <- fit_boot$parameter$theta
      }
      if (save_item_parameters) {
        alpha_store[, , b] <- fit_boot$parameter$alpha
        beta_store[, , b] <- fit_boot$parameter$beta
      }
      if (b %% verbose == 0) {
        cat("* Bootstrap", b, "/", iter, "\n")
      }
      if (verbose <= iter & b == 1) {
        cat("* Bootstrap", 1, "/", iter, "\n")
      }
    }
  } else {
    cl <- parallel::makeCluster(thread)
    doParallel::registerDoParallel(cl)
    if (!is.null(seed)) {
      doRNG::registerDoRNG(seed)
    } else {
      doRNG::registerDoRNG(1)
    }
    fit_foreach <- foreach::foreach(
      b = 1:iter, 
      .packages = "poEMirt"
    ) %dorng% {
      # Simulated response
      dd <- fit$info$data
      dd$data$response <- predict.poEMirtFit(fit, type = "response", return = "array")
      dd$data$modelinput <- construct_sb_auxs(dd$data$response, dd$data$trial, dd$categories)
      fit_boot <- quiet(
        poEMirt(
          data = dd,
          model = fit$info$model,
          constraint = fit$info$constraint,
          alpha_fix = fit$info$alpha_fix,
          theta_std = fit$info$theta_std,
          priors = fit$info$priors,
          control = fit$info$control
        )
      )
      LL <- list()
      LL$theta <- fit_boot$parameter$theta
      if (save_item_parameters) {
        LL$alpha <- fit_boot$parameter$alpha
        LL$beta <- fit_boot$parameter$beta
      }
      return(LL)
    }
    for (b in 1:iter) {
      if (fit$info$model == "static") {
        theta_store[, b] <- fit_foreach[[b]]$theta
      } else {
        theta_store[, , b] <- fit_foreach[[b]]$theta
      }
      if (save_item_parameters) {
        alpha_store[, , b] <- fit_foreach[[b]]$alpha
        beta_store[, , b] <- fit_foreach[[b]]$beta
      }
    }
    parallel::stopCluster(cl)
    rm(list = "fit_foreach")
  }
  L <- list(
    theta = theta_store
  )
  if (save_item_parameters) {
    L$alpha <- alpha_store
    L$beta <- beta_store
  } else {
    L$alpha <- NULL
    L$beta <- NULL
  }
  return(L)
}