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
  sums <- colSums(data$response, na.rm = TRUE)
  for (j in 1:data$size$J) {
    #Kj <- max(data$categories[[j]])
    cat <- data$categories[[j]]
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
  Yimp <- data$response[, , 1] / data$trial
  Yimp <- apply(Yimp, 2, med_impute)
  theta_init <- scale(stats::prcomp(Yimp)$x[, 1])
  if (is.null(constraint)) {
    if (exists('dynamic', data)) {
      constraint <- rep(which.max(theta_init), data$size$T)
    } 
  }
  if (exists('dynamic', data)) {
    # Static IRT
    static_fit <- poEMirtbase_fit(
      Y = data$response,
      N = data$trial,
      alpha_init = alpha_init,
      beta_init = beta_init,
      theta_init = theta_init,
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
  if (fit$info$model == 'static') {
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
      dd$response <- predict.poEMirtFit(fit, type = 'response')
      fit_boot <- quiet(
        poEMirt(
          data = dd,
          model = fit$info$model,
          constraint = fit$info$constraint,
          priors = fit$info$priors,
          control = fit$info$control
        )
      )
      if (fit$info$model == 'static') {
        theta_store[, b] <- fit_boot$parameter$theta
      } else {
        theta_store[, , b] <- fit_boot$parameter$theta
      }
      if (save_item_parameters) {
        alpha_store[, , b] <- fit_boot$parameter$alpha
        beta_store[, , b] <- fit_boot$parameter$beta
      }
      if (b %% verbose == 0) {
        cat('* Bootstrap', b, '/', iter, '\n')
      }
      if (verbose <= iter & b == 1) {
        cat('* Bootstrap', 1, '/', iter, '\n')
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
      .packages = 'poEMirt'
    ) %dorng% {
      # Simulated response
      dd <- fit$info$data
      dd$response <- poEMirt:::predict.poEMirtFit(fit, type = 'response')
      fit_boot <- poEMirt:::quiet(
        poEMirt::poEMirt(
          data = dd,
          model = fit$info$model,
          constraint = fit$info$constraint,
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
      if (fit$info$model == 'static') {
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
    rm(list = 'fit_foreach')
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

#' @description Data-simulation function
#' @importFrom dplyr %>% as_tibble mutate row_number bind_rows tibble left_join relocate arrange
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_remove
#' @importFrom stats qlogis plogis runif rnorm rmultinom 
#' @keywords internal
#' @noRd
SimulateData <- function(I, J, maxK, T = NULL, n_range) {
  if (!is.null(T)) {
    # Timemaps
    timemap <- matrix(1, I, T)
    item_timemap <- rep(0:(T-1), each = J/T)
    item_match <- rep(NA, J)
  }
  # Number of choice
  if(maxK == 2) {
    K <- rep(maxK, J)
  } else {
    K <- sample(2:maxK, J, replace = TRUE)
  }
  # Item parameters
  alpha <- beta <- matrix(NA, J, maxK - 1)
  for (j in 1:J) {
    Kj <- K[j]
    baseline_probs <- runif(Kj)
    baseline_probs <- exp(baseline_probs)/sum(exp(baseline_probs))
    baseline_sb <- cumprod(1-baseline_probs)/(1-baseline_probs) * baseline_probs
    baseline_sb <- baseline_sb[1:(Kj-1)] 
    baseline_sb <- stats::qlogis(baseline_sb)
    x1_probs <- stats::runif(Kj)
    x1_probs <- exp(x1_probs)/sum(exp(x1_probs))
    x1_sb <- cumprod(1-x1_probs)/(1-x1_probs) * x1_probs
    x1_sb <- x1_sb[1:(Kj-1)]    
    x1_sb <- stats::qlogis(x1_sb)
    alpha[j, 1:(Kj-1)] <- baseline_sb
    beta[j, 1:(Kj-1)] <- x1_sb - baseline_sb
  }
  if (!is.null(T)) {
    # Latent traits
    theta <- matrix(NA, I, T)
    theta[, 1] <- stats::rnorm(I)
    for (t in 2:T) theta[, t] <- stats::rnorm(I, theta[, t-1], 0.1)
    con <- apply(theta, 2, which.max)
  } else {
    theta <- stats::rnorm(I)
    con <- which.max(theta)
  }
  
  ## Responses
  if (length(n_range) == 1) {
    n <- matrix(n_range, I, J)
  } else {
    n <- sample(n_range[1]:n_range[2], size = I * J, replace = TRUE) %>% 
      matrix(I, J)
  }
  Y <- array(NA, c(I, J, maxK))
  for (i in 1:I) {
    for (j in 1:J) {
      Kj <- K[j]
      if (!is.null(T)) t <- item_timemap[j] + 1
      if (Kj > 2) {
        if (!is.null(T)) {
          psi <- alpha[j, 1:(Kj-1)] + beta[j, 1:(Kj-1)] * theta[i, t]
        } else {
          psi <- alpha[j, 1:(Kj-1)] + beta[j, 1:(Kj-1)] * theta[i]
        }
        psi <- stats::plogis(psi)
        p <- psi * c(1, cumprod(1 - psi))[-length(psi)]
        p <- c(p, 1 - sum(p))
      } else {
        if (!is.null(T)) {
          psi <- alpha[j, 1] + beta[j, 1] * theta[i, t]
        } else {
          psi <- alpha[j, 1] + beta[j, 1] * theta[i]
        }
        psi <- stats::plogis(psi)
        p <- c(psi, 1 - psi)
      }
      Y[i, j, 1:Kj] <- stats::rmultinom(1, n[i, j], p)
    }
  }
  
  # Check
  categories <- apply(colSums(Y, na.rm = TRUE), 1, function(x) which(x != 0), simplify = FALSE)
  check <- lapply(
    categories,
    function(x) {
      all(sort(x) == order(sort(x)))
    }
  ) %>% 
    unlist() %>% 
    sum()
  if (check != J) message('* WARNING: Some items have fewer response categories than your setting.')
  
  # Converting into dataframe 
  d <- rep()
  for (k in 1:maxK) {
    if (k == 1) {
      d <- Y[, , k] %>% 
        dplyr::as_tibble() %>% 
        suppressWarnings() %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -.data$i,
          names_to = 'j',
          values_to = paste0('y', k)
        ) %>% 
        dplyr::mutate(
          j = .data$j %>% 
            stringr::str_remove('V') %>% 
            as.integer()
        ) %>% 
        dplyr::bind_rows(d, .)
    } else {
      d <- Y[, , k] %>% 
        dplyr::as_tibble() %>% 
        suppressWarnings() %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -.data$i,
          names_to = 'j',
          values_to = paste0('y', k)
        ) %>% 
        dplyr::select(-i, -j) %>% 
        dplyr::bind_cols(d, .)
    }
  }
  if (!is.null(T)) {
    d <- d %>% 
      dplyr::left_join(
        dplyr::tibble(
          j = 1:J,
          t = item_timemap + 1
        ),
        by = 'j'
      ) %>% 
      dplyr::relocate(.data$t, .before = .data$y1) %>% 
      dplyr::arrange(.data$j)
  }
  
  L <- list(
    df = d,
    parameter = list(
      theta = theta,
      alpha = alpha,
      beta = beta
    )
  )
  return(L)
}


