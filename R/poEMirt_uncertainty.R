#' @title Estimating statistical uncertainty for poEMirt models
#' @description This function computes statistical uncertainty using bootstrap or Gibbs sampling.
#' @param fit A \code{poEMirtFit} object from \code{poEMirt()}.
#' @param method A character. This must be "bootstrap" or "gibbs". 
#' @param seed An integer. Random seed if needed.
#' @param iter An integer. The number of iterations. Default is 100.
#' @param control A list of model controls. This can include following elements.
#' \itemize{
#'   \item \code{save_item_parameters} A bool. If TRUE, the function keeps item parameters. Default is TRUE.
#'   \item \code{verbose} An integer of verbose. Default is NULL.
#'   \item \code{thread} An integer of threads (only for Bootstrap). Default is 1. Status messages will not be printed if \code{thread} > 1.
#'   \item \code{warmup} An integer (only for Gibbs). The number of warmup iterations for Gibbs sampling. Default is 100.
#'   \item \code{thin} An integer of thinning for Gibbs sampling (only for Gibbs). Default is 1.
#'   \item \code{PG_approx} A bool (only Gibbs). If TRUE, approximated Polya-Gamma random draws are used. This approximation accelerates the implementation but is not recommended for small n (e.g., n < 10). Default is FALSE.
#' }
#' @return A \code{poEMirtBoot} / \code{poEMirtGibbs} object containing:
#' \describe{
#'   \item{parameter}{A list of parameters}
#'   \item{standard_deviation}{A list of standard deviations of parameters}
#'   \item{input}{A list of model input}
#' }
#' @importFrom Rcpp sourceCpp
#' @importFrom stats sd
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
#' # Bootstrap
#' fit_boot <- poEMirt_uncertainty(
#'   fit = fit,
#'   method = "bootstrap",
#'   seed = 1,
#'   iter = 100,
#'   control = list(
#'     verbose = 10
#'   )
#' )
#' summary(fit_boot, parameter = "theta", ci = 0.95)
#' 
#' # Gibbs
#' fit_gibbs <- poEMirt_uncertainty(
#'   fit = fit,
#'   method = "gibbs",
#'   seed = 1,
#'   iter = 500,
#'   control = list(
#'     verbose = 50,
#'     PG_approx = TRUE, 
#'     warmup = 100,
#'     thin = 5
#'   )
#' )
#' summary(fit_gibbs, parameter = "theta", ci = 0.95)
#' }

poEMirt_uncertainty <- function(fit, 
                                method = c("bootstrap", "gibbs"), 
                                seed = NULL,
                                iter = 100, 
                                control = NULL) 
{
  if (class(fit)[1] != "poEMirtFit") stop("`fit` must be 'poEMirtFit' object from poEMirt::poEMirt().")
  method <- match.arg(method, choices = c("bootstrap", "gibbs"))
  if (!is.null(seed)) set.seed(seed)
  if (method == "bootstrap") {
    # Bootstrap
    cat("=== Parametric bootstrap to estimate statistical uncertainty for poEMirt ===\n")
    if (is.null(control)) control <- list()
    if (!exists("save_item_parameters", control)) control$save_item_parameters <- TRUE
    if (!exists("verbose", control)) control$verbose <- NULL
    if (!exists("thread", control)) {
      control$thread <- 1
    } else {
      if (control$thread != 1) {
        cat("* Computing with", control$thread, "threads.....\n")
      }
    }
    verb <- ifelse(is.null(control$verbose), iter+1, control$verbose)
    stime <- proc.time()[3]
    fit_2 <- poEMirt_boot(
      fit = fit,
      iter = iter,
      verbose = verb,
      save_item_parameters = control$save_item_parameters,
      thread = control$thread,
      seed = seed
    )
    etime <- proc.time()[3]
    el <- round(etime - stime, 1)
  } else {
    # Gibbs sampling
    cat("=== Gibbs Sampling to estimate statistical uncertainty for poEMirt ===\n")
    if (is.null(control)) control <- list()
    if (!exists("PG_approx", control)) control$PG_approx <- FALSE
    if (control$PG_approx) {
      meann <- mean(fit$info$data$data$trial)
      if (meann < 10) {
        message("* NOTE: `PG_approx = TRUE` is not recommended for small number of trials (e.g., n < 10)")
      }
    }
    if (!exists("warmup", control)) control$warmup <- 100
    if (!exists("thin", control)) control$thin <- 1
    if (!exists("save_item_parameters", control)) control$save_item_parameters <- TRUE
    if (!exists("verbose", control)) control$verbose <- NULL
    verb <- ifelse(is.null(control$verbose), iter+1, control$verbose)
    fit$parameter$alpha[is.na(fit$parameter$alpha)] <- 0
    fit$parameter$beta[is.na(fit$parameter$beta)] <- 0
    fit$parameter$theta[is.na(fit$parameter$theta)] <- 0
    stime <- proc.time()[3]
    if (fit$info$model == "static") {
      if (is.null(control$verbose)) {
        fit_2 <- quiet(
          poEMirtbase_gibbs_fit(
            Y = fit$info$data$data$response, 
            S = fit$info$data$data$modelinput$S, 
            Nks = fit$info$data$data$modelinput$Nks, 
            alpha = fit$parameter$alpha, 
            beta = fit$parameter$beta, 
            theta = fit$parameter$theta, 
            unique_categories = fit$info$data$categories, 
            a0 = fit$info$priors$a0,
            A0 = fit$info$priors$A0,
            b0 = fit$info$priors$b0,
            B0 = fit$info$priors$B0,
            PG_approx = control$PG_approx,
            constraint = fit$info$constraint - 1,
            std = fit$info$theta_std,
            iter = iter,
            warmup = control$warmup,
            thin = control$thin,
            save_item_parameters = control$save_item_parameters,
            verbose = verb
          )
        )
      } else {
        fit_2 <- 
          poEMirtbase_gibbs_fit(
            Y = fit$info$data$data$response, 
            S = fit$info$data$data$modelinput$S, 
            Nks = fit$info$data$data$modelinput$Nks, 
            alpha = fit$parameter$alpha, 
            beta = fit$parameter$beta, 
            theta = fit$parameter$theta, 
            unique_categories = fit$info$data$categories, 
            a0 = fit$info$priors$a0,
            A0 = fit$info$priors$A0,
            b0 = fit$info$priors$b0,
            B0 = fit$info$priors$B0,
            PG_approx = control$PG_approx,
            constraint = fit$info$constraint - 1,
            std = fit$info$theta_std,
            iter = iter,
            warmup = control$warmup,
            thin = control$thin,
            save_item_parameters = control$save_item_parameters,
            verbose = verb
          )
      }
    } else {
      # Dynamic
      if (is.null(control$verbose)) {
        fit_2 <- quiet(
          poEMirtdynamic_gibbs_fit(
            Y = fit$info$data$data$response, 
            S = fit$info$data$data$modelinput$S, 
            Nks = fit$info$data$data$modelinput$Nks, 
            sb_check = fit$info$data$data$modelinput$sb_check, 
            alpha = fit$parameter$alpha, 
            beta = fit$parameter$beta, 
            theta = fit$parameter$theta, 
            unique_categories = fit$info$data$categories, 
            uJ_J = fit$info$data$uniqueJ_J, 
            timemap2 = fit$info$data$dynamic$timemap2, 
            item_timemap = fit$info$data$dynamic$item_timemap, 
            IT = fit$info$data$dynamic$index$IT, 
            ITJ = fit$info$data$dynamic$index$ITJ, a0 = fit$info$priors$a0,
            A0 = fit$info$priors$A0,
            b0 = fit$info$priors$b0,
            B0 = fit$info$priors$B0,
            m0 = fit$info$priors$m0,
            C0 = fit$info$priors$C0,
            Delta = fit$info$priors$Delta,
            alpha_fix = fit$info$alpha_fix, 
            PG_approx = control$PG_approx, 
            constraint = fit$info$constraint - 1, 
            std = fit$info$theta_std, 
            iter = iter,
            warmup = control$warmup,
            thin = control$thin,
            save_item_parameters = control$save_item_parameters,
            verbose = verb
          )
        )
      } else {
        fit_2 <- poEMirtdynamic_gibbs_fit(
          Y = fit$info$data$data$response, 
          S = fit$info$data$data$modelinput$S, 
          Nks = fit$info$data$data$modelinput$Nks, 
          sb_check = fit$info$data$data$modelinput$sb_check, 
          alpha = fit$parameter$alpha, 
          beta = fit$parameter$beta, 
          theta = fit$parameter$theta, 
          unique_categories = fit$info$data$categories, 
          uJ_J = fit$info$data$uniqueJ_J, 
          timemap2 = fit$info$data$dynamic$timemap2, 
          item_timemap = fit$info$data$dynamic$item_timemap, 
          IT = fit$info$data$dynamic$index$IT, 
          ITJ = fit$info$data$dynamic$index$ITJ, 
          a0 = fit$info$priors$a0,
          A0 = fit$info$priors$A0,
          b0 = fit$info$priors$b0,
          B0 = fit$info$priors$B0,
          m0 = fit$info$priors$m0,
          C0 = fit$info$priors$C0,
          Delta = fit$info$priors$Delta,
          alpha_fix = fit$info$alpha_fix, 
          PG_approx = control$PG_approx, 
          constraint = fit$info$constraint - 1, 
          std = fit$info$theta_std, 
          iter = iter,
          warmup = control$warmup,
          thin = control$thin,
          save_item_parameters = control$save_item_parameters,
          verbose = verb
        )
      }
    }
  }
  etime <- proc.time()[3]
  el <- round(etime - stime, 1)
  cat("* DONE!\n")
  sds <- list()
  if (fit$info$model == "dynamic") {
    if (method == "gibbs") {
      th <- array(NA, c(fit$info$data$size$I, fit$info$data$size$T, iter/control$thin))
      rownames(th) <- rownames(fit$parameter$theta)
      colnames(th) <- colnames(fit$parameter$theta)
      for (ii in 1:(iter/control$thin)) {
        th[, , ii] <- fit_2$theta[[ii]]
      }
      th[th == 0] <- NA
      fit_2$theta <- th
      rm(list = "th")
    } else {
      dimnames(fit_2$theta) <- list(rownames(fit$parameter$theta), colnames(fit$parameter$theta))
    }
    sds$theta <- apply(fit_2$theta, c(1, 2), stats::sd, na.rm = TRUE)
    dimnames(sds$theta) <- list(rownames(fit$parameter$theta), colnames(fit$parameter$theta))
  } else {
    rownames(fit_2$theta) <- rownames(fit$parameter$theta)
    colnames(fit_2$theta) <- colnames(fit$parameter$theta)
    sds$theta <- apply(fit_2$theta, 1, stats::sd, na.rm = TRUE)
    fit_2$theta[fit_2$theta == 0] <- NA
    dimnames(sds$theta) <- list(rownames(fit$parameter$theta), colnames(fit$parameter$theta))
  }
  
  if (control$save_item_parameters) {
    if (method == "gibbs") {
      al <- be <- array(NA, c(fit$info$data$size$J, fit$info$data$size$maxK-1, iter/control$thin))
      rownames(al) <- rownames(be) <- rownames(fit$parameter$alpha)
      colnames(al) <- colnames(be) <- colnames(fit$parameter$alpha)
      for (ii in 1:(iter/control$thin)) {
        al[, , ii] <- fit_2$alpha[[ii]]
        be[, , ii] <- fit_2$beta[[ii]]
      }
      al[al == 0] <- NA
      be[be == 0] <- NA
      fit_2$alpha <- al
      fit_2$beta <- be
      rm(list = c("al", "be"))
    } else {
      rownames(fit_2$alpha) <- rownames(fit_2$beta) <- rownames(fit$parameter$alpha)
      colnames(fit_2$alpha) <- colnames(fit_2$beta) <- colnames(fit$parameter$beta)
    }
    sds$alpha <- apply(fit_2$alpha, c(1, 2), stats::sd, na.rm = TRUE)
    sds$beta <- apply(fit_2$beta, c(1, 2), stats::sd, na.rm = TRUE)
    dimnames(sds$alpha) <- dimnames(sds$beta) <- list(rownames(fit$parameter$alpha), colnames(fit$parameter$beta))
  }
  L <- list(
    parameter = fit_2,
    standard_deviation = sds,
    input = list(
      fit = fit,
      method = method,
      iter = iter,
      control = control
    )
  )
  L$call <- match.call()
  L$time <- list(
    date = date(),
    start = stime,
    end = etime,
    elapsed = el
  )
  class(L) <- c(paste0("poEMirt", ifelse(method == "bootstrap", "Boot", "Gibbs")), class(L))
  
  return(L)
}