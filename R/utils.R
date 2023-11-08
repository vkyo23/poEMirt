#' @description Imputing function
#' @importFrom stats median
#' @keywords internal
#' @noRd
med_impute <- function(x) {
  x[is.na(x) | is.infinite(x)] <- stats::median(x, na.rm = TRUE)
  x[is.na(x) | is.infinite(x)] <- 0
  return(x)
}

#' @keywords internal
#' @useDynLib poEMirt, .registration = TRUE
#' @noRd
predict1 <- function(object,
                     return_array = FALSE,
                     type = "prob") 
{
  object$parameter$alpha[is.na(object$parameter$alpha)] <- 0
  object$parameter$beta[is.na(object$parameter$beta)] <- 0
  object$parameter$theta[is.na(object$parameter$theta)] <- 0
  if (object$info$model != "static") {
    out <- prediction(
      Y = object$info$data$data$response,
      N = object$info$data$data$trial,
      alpha = object$parameter$alpha,
      beta = object$parameter$beta,
      theta = object$parameter$theta,
      unique_categories = object$info$data$categories,
      item_timemap = object$info$data$dynamic$item_timemap,
      model = "dynamic",
      type = type
    )
  } else {
    out <- prediction(
      Y = object$info$data$data$response,
      N = object$info$data$data$trial,
      alpha = object$parameter$alpha,
      beta = object$parameter$beta,
      theta = object$parameter$theta,
      unique_categories = object$info$data$categories,
      item_timemap = NA,
      model = "static",
      type = type
    )
  }
  dimnames(out) <- dimnames(object$info$data$data$response)
  if (!return_array) {
    maxK <- dim(out)[3]
    for (k in 1:maxK) {
      nam <- dimnames(out)[[3]][k]
      tmp1 <- out[, , nam] %>%
        dplyr::as_tibble(.name_repair = "unique") %>%
        dplyr::mutate(i = rownames(out)) %>%
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "j",
          values_to = nam
        )
      if (k == 1) {
        tmp <- tmp1
      } else {
        tmp <- dplyr::bind_cols(tmp, tmp1[, 3])
      }
    }
    tmp <- tmp %>%
      dplyr::mutate(type, .before = "j")
    incl <- which(apply(tmp[, -c(1:3)], 1, function(x) sum(is.na(x))) != ncol(tmp[, -c(1:3)]))
    tmp <- tmp[incl, ]
    if (exists("rep", object$info$data)) {
      out <- tmp %>%
        dplyr::mutate(
          index = stringr::str_c("[", .data$i, ",", stringr::str_replace(.data$j, "-", ","),"]"),
          reference = "[i,t,j]"
        ) %>%
        dplyr::select(-"i", -"j") %>%
        dplyr::relocate("index", "reference", .after = "type")
    } else {
      out <- tmp %>%
        dplyr::mutate(
          index = stringr::str_c("[", .data$i, ",", .data$j,"]"),
          reference = "[i,j]"
        ) %>%
        dplyr::select(-"i", -"j") %>%
        dplyr::relocate("index", "reference", .after = "type")
    }
  }
  
  return(out)
}

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
    if (fit$info$estimate_Delta) {
      Delta_store <- rep()
    }
  }
  if (save_item_parameters) {
    alpha_store <- beta_store <- array(NA, c(fit$info$data$size$J, fit$info$data$size$maxK-1, iter))
  }
  if (thread == 1) {
    for (b in 1:iter) {
      # Simulated response
      dd <- fit$info$data
      dd$data$response <- predict.poEMirtFit(fit, return_array = TRUE, type = "response")
      dd$data$modelinput <- construct_sb_auxs(dd$data$response, dd$data$trial, dd$categories)
      fit_boot <- quiet(
        poEMirt(
          data = dd,
          model = fit$info$model,
          constraint = fit$info$constraint,
          fix_alpha = fit$info$fix_alpha,
          fix_beta = fit$info$fix_beta,
          estimate_Delta = fit$info$estimate_Delta,
          std_theta = fit$info$std_theta,
          priors = fit$info$priors,
          control = fit$info$control
        )
      )
      if (fit$info$model == "static") {
        theta_store[, b] <- fit_boot$parameter$theta
      } else {
        theta_store[, , b] <- fit_boot$parameter$theta
        if (fit$info$estimate_Delta) {
          Delta_store[b] <- fit_boot$parameter$Delta
        }
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
      dd$data$response <- predict.poEMirtFit(fit, return_array = TRUE, type = "response")
      dd$data$modelinput <- construct_sb_auxs(dd$data$response, dd$data$trial, dd$categories)
      fit_boot <- quiet(
        poEMirt(
          data = dd,
          model = fit$info$model,
          constraint = fit$info$constraint,
          fix_alpha = fit$info$fix_alpha,
          fix_beta = fit$info$fix_beta,
          estimate_Delta = fit$info$estimate_Delta,
          std_theta = fit$info$std_theta,
          priors = fit$info$priors,
          control = fit$info$control
        )
      )
      LL <- list()
      LL$theta <- fit_boot$parameter$theta
      if (fit$info$estimate_Delta) {
        LL$Delta <- fit_boot$parameter$Delta
      }
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
        if (fit$info$estimate_Delta) {
          Delta_store[b] <- fit_foreach[[b]]$Delta
        }
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
  if (fit$info$model == "dynamic" & fit$info$estimate_Delta) {
    L$Delta <- Delta_store
  }
  if (save_item_parameters) {
    L$alpha <- alpha_store
    L$beta <- beta_store
  } else {
    L$alpha <- NULL
    L$beta <- NULL
  }
  return(L)
}

