#' @title A fast EM Item Response Theory model for public opinion analysis
#' @description Fit \code{poEMirt} models.
#' 
#' @param data An object of \code{poEMirtData} via \code{read_poEMirt()}. 
#' @param model A string, one of "static" or "dynamic".
#' @param constraint An integer scalar or vector (for dynamic model, the same length as the number of times), index of an individual i (the location of i) whose latent trait is always set positive.
#' @param alpha_fix A bool, whether fixes alpha of same repeated items or not. Default is FALSE.
#' @param theta_std A bool, whether standardizes theta or not. Default is FALSE
#' @param init A list, containing initial values (optional).
#' \itemize{
#'   \item \code{alpha} J x max(K)-1 matrix of alpha. Missing values must be 0.
#'   \item \code{beta} J x max(K)-1 matrix of beta. Missing values must be 0.
#'   \item \code{theta} I-length vector (for static model) or I x T matrix of theta (for dynamic model). Missing values must be 0.
#' }
#' @param priors a list, containing prior distributions (optional).
#' \itemize{
#'   \item \code{a0} A double or J x max(K)-1 matrix, prior mean of alpha. Default is 0.
#'   \item \code{A0} A double or J x max(K)-1 matrix, prior variance of alpha. Default is 1.
#'   \item \code{b0} A double or J x max(K)-1 matrix, prior mean of beta Default is 0.
#'   \item \code{B0} A double or J x max(K)-1 matrix, prior variance of beta. Default is 1.
#'   \item \code{m0} A double or I length vector, prior mean of theta_i0 for dynamic model. Default is 0. 
#'   \item \code{C0} A double or I length vector, prior variance of theta_i0 for dynamic model. Default is 1.
#'   \item \code{Delta} A double or I length vector, prior evolution variance of theta_it for dynamic model. Default is 0.01.
#' }
#' @param control A list of model controls.
#' \itemize{
#'   \item \code{compute_ll} A bool, whether compute log-likelihood for each iteration or not. Default is FALSE.
#'   \item \code{maxit} An integer, maximum number of iterations for EM. Default is 500.
#'   \item \code{tol} A double (< 1), convergence threshold. Default is 1e-6.
#'   \item \code{verbose} An integer, the function prints the status every `verbose`.
#' }
#' @return A \code{poEMirtFit} object containing:
#' \describe{
#'   \item{parameter}{A list of parameters}
#'   \item{iteration}{The number of iterations (integer)}
#'   \item{converge}{Whether the model converged or not (boolean)}
#'   \item{convstat}{A matrix of convergence statistics}
#'   \item{info}{A list of model information}
#' }
#' @importFrom Rcpp sourceCpp
#' @rdname poEMirt
#' 
#' @useDynLib poEMirt, .registration = TRUE
#' @export
#' 
#' @examples 
#' \dontrun{
#' data("sim_data_dynamic")
#' 
#' # Convert into poEMirt-readable data
#' data <- read_poEMirt(
#'   data = sim_data_dynamic,
#'   responses = paste0('y', 1:5),
#'   i = "i",
#'   j = "j",
#'   t = "t"
#' )
#' 
#' # Fit the model
#' fit <- poEMirt(
#'   data = data,
#'   model = "dynamic",
#'   control = list(
#'     verbose = 10,
#'     constrant = 1
#'   )
#' )
#' 
#' # Summarize the result 
#' summary(fit, parameter = "theta")
#' }

poEMirt <- function(data, 
                    model = c("static", "dynamic"), 
                    constraint = NULL,
                    alpha_fix = FALSE,
                    theta_std = FALSE,
                    init = NULL, 
                    priors = NULL, 
                    control = NULL) 
{
  cat("=== poEMirt starts! ===\n")
  
  # Input check
  if (is.null(control)) control <- list()
  model <- match.arg(model, choices = c("static", "dynamic"))
  if (length(model) > 1) stop("Model must be one of 'static' or 'dynamic'.")
  if (class(data)[1] != "poEMirtData") {
    stop("`data` should be a `poEMirtData` object. Use `read_poEMirt()` first to create the object.")
  }
  I <- data$size$I
  J <- data$size$J
  maxK <- data$size$maxK
  if (model == "dynamic") {
    T <- data$size$T
    timemap <- data$dynamic$timemap
    item_timemap <- data$dynamic$item_timemap
  }
  unq_cat <- data$categories
  
  if (alpha_fix) {
    if (!exists("rep", data)) {
      stop("Cannot detect repeated items. Try `alpha_fix = FALSE` or make sure your j in `read_poEMirt()`.")
    } else {
      uJ_J <- data$rep$processed
    }
  } else {
    uJ_J <- list(NA)
  }
  
  # Priors
  if (is.null(priors)) {
    cat("* Setting priors.....")
    priors <- list(
      a0 = matrix(0, J, maxK-1),
      A0 = matrix(1, J, maxK-1),
      b0 = matrix(0, J, maxK-1),
      B0 = matrix(1, J, maxK-1)
    )
    if (model == "dynamic") {
      priors$m0 <- rep(0, I)
      priors$C0 <- rep(1, I)
      priors$Delta <- rep(0.01, I)
    }
    cat("DONE!\n")
  } else {
    if (exists("a0", priors)) {
      if (length(priors$a0) == 1) {
        priors$a0 <- matrix(priors$a0, J, maxK-1)
      } else if (length(priors$a0) != (J * (maxK-1))) {
        stop("`priors$a0` must be must be a scalar or J x max(Kj)-1 matrix.")
      }
    } else {
      priors$a0 <- matrix(0, J, maxK-1)
    }
    if (exists("A0", priors)) {
      if (length(priors$A0) == 1) {
        priors$A0 <- matrix(priors$A0, J, maxK-1)
      } else if (length(priors$A0) != (J * (maxK-1))) {
        stop("`priors$A0` must be must be a scalar or J x max(Kj)-1 matrix.")
      }
    } else {
      priors$A0 <- matrix(1, J, maxK-1)
    }
    if (exists("b0", priors)) {
      if (length(priors$b0) == 1) {
        priors$b0 <- matrix(priors$b0, J, maxK-1)
      } else if (length(priors$b0) != (J * (maxK-1))) {
        stop("`priors$b0` must be must be a scalar or J x max(Kj)-1 matrix.")
      }
    } else {
      priors$b0 <- matrix(0, J, maxK-1)
    }
    if (exists("B0", priors)) {
      if (length(priors$B0) == 1) {
        priors$B0 <- matrix(priors$B0, J, maxK-1)
      } else if (length(priors$B0) != (J * (maxK-1))) {
        stop("`priors$B0` must be must be a scalar or J x max(Kj)-1 matrix.")
      }
    } else {
      priors$B0 <- matrix(1, J, maxK-1)
    }
    if (model == "dynamic") {
      if (exists("m0", priors)) {
        if (length(priors$m0) == 1) {
          priors$m0 <- rep(priors$m0, I)
        } else if (length(priors$m0) != I) {
          stop("`priors$m0` must be must be a scalar or I-length vector.")
        }
      } else {
        priors$m0 <- rep(0, I)
      }
      if (exists("C0", priors)) {
        if (length(priors$C0) == 1) {
          priors$C0 <- rep(priors$C0, I)
        } else if (length(priors$C0) != I) {
          stop("`priors$C0` must be must be a scalar or I-length vector.")
        }
      } else {
        priors$C0 <- rep(1, I)
      }
      if (exists("Delta", priors)) {
        if (length(priors$Delta) == 1) {
          priors$Delta <- rep(priors$Delta, I)
        } else if (length(priors$Delta) != I) {
          stop("`priors$Delta` must be must be a scalar or I-length vector.")
        }
      } else {
        priors$Delta <- rep(0.01, I)
      }
    } 
  }
  
  # Initial value
  if (is.null(init)) {
    cat("* Finding best initial values.....")
    init <- make_init(
      data = data,
      priors = priors,
      constraint = constraint
    )
    if (model == "dynamic") {
      init$theta <- init$theta * timemap
    }
    cat("DONE!\n")
  } else {
    inls <- c("alpha", "beta", "theta")
    for (l in 1:length(inls)) {
      if (!exists(inls[l], init)) {
        stop(paste0("A list `init` does not contain `", inls[l], "`."))
      } else if ((inls[l] %in% c("alpha", "beta") & nrow((init[inls[l]])[[1]]) != J) |
                 (inls[l] %in% c("alpha", "beta") & ncol(init[inls[l]][[1]]) != (maxK-1))) {
        stop(paste0("`init$", inls[l], "` must be J x K-1 matrix"))
      } else if (inls[l] == "theta" & model == "dynamic") {
        if (nrow(init[inls[l]][[1]]) != I | ncol(init[inls[l]][[1]]) != T) {
          stop("`init$theta` must be I x T matrix for dynamic model.") 
        }
      } else if ((inls[l] == "theta" & model == "static" & nrow(init[inls[l]][[1]]) != I) |
                 (inls[l] == "theta" & model == "static" & ncol(init[inls[l]][[1]]) != 1)) {
        stop("`init$theta` must be I-length vector for static model.")
      }
    }
  }
  # Constraint
  if (is.null(constraint)) {
    if (model == "dynamic") {
      constraint <- apply(init$theta, 2, which.max)
    } else {
      constraint <- which.max(init$theta) 
    }
  } 
  
  # Control
  if (!exists("compute_ll", control)) control$compute_ll <- FALSE
  if (!exists("maxit", control)) control$maxit <- 500
  if (!exists("tol", control)) control$tol <- 1e-6
  if (control$tol >= 1) stop("`control$tol` must be lower than 1.")
  if (!exists("verbose", control)) control$verbose <- NULL
  verb <- ifelse(is.null(control$verbose), control$maxit+1, control$verbose)
  
  # Fitting
  cat("* Expectation-Maximization\n")
  stime <- proc.time()[3]
  if (model == "dynamic") {
    fit <- poEMirtdynamic_fit(
      Y = data$data$response, 
      S = data$data$modelinput$S, 
      Nks = data$data$modelinput$Nks, 
      alpha_old = init$alpha, 
      beta_old = init$beta, 
      theta_old = init$theta, 
      unique_categories = data$categories, 
      uJ_J = uJ_J, 
      timemap2 = data$dynamic$timemap2, 
      item_timemap = item_timemap, 
      IT = data$dynamic$index$IT, 
      ITJ = data$dynamic$index$ITJ, 
      a0 = priors$a0, 
      A0 = priors$A0, 
      b0 = priors$b0, 
      B0 = priors$B0, 
      m0 = priors$m0, 
      C0 = priors$C0, 
      Delta = priors$Delta, 
      constraint = constraint - 1, 
      alpha_fix = alpha_fix, 
      std = theta_std, 
      maxit = control$maxit, 
      verbose = verb, 
      tol = control$tol, 
      compute_ll = control$compute_ll
    )
  } else {
    fit <- poEMirtbase_fit(
      Y = data$data$response, 
      S = data$data$modelinput$S, 
      Nks = data$data$modelinput$Nks, 
      alpha_old = init$alpha, 
      beta_old = init$beta, 
      theta_old = init$theta, 
      unique_categories = data$categories, 
      a0 = priors$a0, 
      A0 = priors$A0, 
      b0 = priors$b0, 
      B0 = priors$B0, 
      constraint = constraint - 1, 
      std = theta_std, 
      maxit = control$maxit, 
      verbose = verb, 
      tol = control$tol, 
      compute_ll = control$compute_ll
    )
  }
  # Output
  etime <- proc.time()[3]
  el <- round(etime - stime, 1)
  if (fit$converge) {
    cat("* Model converged at iteration", fit$iter, ":", el, "sec.\n")
  } else {
    warning("* Model failed to converge : ", el, "sec.")
  }
  fit$alpha[fit$alpha == 0] <- NA
  fit$beta[fit$beta == 0] <- NA
  fit$theta[fit$theta == 0] <- NA
  rownames(fit$alpha) <- rownames(fit$beta) <- colnames(data$data$response)
  colnames(fit$alpha) <- colnames(fit$beta) <- dimnames(data$data$response)[[3]][-maxK]
  rownames(fit$theta) <- rownames(data$data$response)
  if (model == "dynamic") {
    colnames(fit$theta) <- 1:data$size$T
  } else {
    colnames(fit$theta) <- 1
  }
  colnames(fit$conv) <- c("alpha", "beta", "theta")
  data$uniqueJ_J <- uJ_J
  L <- list(
    parameter = list(
      alpha = fit$alpha,
      beta = fit$beta,
      theta = fit$theta
    ),
    iteration = fit$iter,
    converge = fit$converge,
    convstat = fit$conv,
    log_likelihood = fit$log_likelihood,
    info = list(
      data = data, 
      model = model, 
      constraint = constraint,
      alpha_fix = alpha_fix,
      theta_std = theta_std,
      init = init, 
      priors = priors, 
      control = control
    )
  )
  if (!control$compute_ll) L$log_likelihood <- NULL
  L$call <- match.call()
  L$time <- list(
    date = date(),
    start = stime,
    end = etime,
    elapsed = el
  )
  class(L) <- c("poEMirtFit", class(L))
  return(L)
}


