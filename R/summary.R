#' @title Summary function for \code{poEMirtFit}
#' @description \code{summary.poEMirtFit} returns a dataframe of estimated parameters.
#' @param object An object class of \code{poEMirtFit}.
#' @param parameter A character or character vector, indicating what parameters ("alpha", "beta", "theta") you want to get. 
#' @param ... Parameters to \code{summary()}.
#' @returns A tibble dataframe.
#' @importFrom dplyr %>% mutate select as_tibble bind_rows row_number
#' @importFrom stringr str_remove str_c
#' @importFrom tidyr pivot_longer drop_na
#' @export
summary.poEMirtFit <- function(object, 
                               parameter = c('alpha', 'beta', 'theta'),
                               ...) {
  out <- rep()
  message('* Summarizing following parameters: ', paste(parameter, collapse = ', '), '.')
  if (any(parameter == 'alpha')) {
    out <- object$parameter$alpha %>% 
      dplyr::as_tibble() %>% 
      dplyr::mutate(j = dplyr::row_number()) %>% 
      tidyr::pivot_longer(
        cols = -j,
        names_to = 'k',
        values_to = 'alpha'
      ) %>% 
      tidyr::drop_na() %>% 
      dplyr::mutate(
        k = k %>% 
          stringr::str_remove('y') %>% 
          as.integer(),
        index = stringr::str_c('[', j, ', ', k, ']'),
        reference = '[j, k]',
        parameter = 'alpha'
      ) %>% 
      dplyr::select(parameter, index, reference, estimate = alpha) %>% 
      dplyr::bind_rows(out, .)
  }
  if (any(parameter == 'beta')) {
    out <- object$parameter$beta %>% 
      dplyr::as_tibble() %>% 
      dplyr::mutate(j = dplyr::row_number()) %>% 
      tidyr::pivot_longer(
        cols = -j,
        names_to = 'k',
        values_to = 'beta'
      ) %>% 
      tidyr::drop_na() %>% 
      dplyr::mutate(
        k = k %>% 
          stringr::str_remove('y') %>% 
          as.integer(),
        index = stringr::str_c('[', j, ', ', k, ']'),
        reference = '[j, k]',
        parameter = 'beta'
      ) %>% 
      dplyr::select(parameter, index, reference, estimate = beta) %>% 
      dplyr::bind_rows(out, .)
  }
  if (any(parameter == 'theta')) {
    tmp <- object$parameter$theta %>% 
      dplyr::as_tibble()
    if (object$info$model == 'static') {
      out <- tmp %>% 
        dplyr::mutate(index = stringr::str_c('[', dplyr::row_number(), ']')) %>% 
        dplyr::mutate(
          reference = '[i]',
          parameter = 'theta'
        ) %>% 
        dplyr::select(parameter, index, reference, estimate = V1) %>% 
        dplyr::bind_rows(out, .)
    } else {
      out <- tmp %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -i,
          names_to = 't',
          values_to = 'theta'
        ) %>% 
        dplyr::mutate(
          t = t %>% 
            stringr::str_remove('V') %>% 
            as.integer(),
          index = stringr::str_c('[', i, ', ', t, ']'),
          reference = '[i, t]',
          parameter = 'theta'
        ) %>% 
        dplyr::select(parameter, index, reference, estimate = theta) %>% 
        dplyr::bind_rows(out, .)
    }
  }
  return(out)
}

#' @title Summary function for \code{poEMirtBoot}
#' @description \code{summary.poEMirtBoot} returns a dataframe of estimated parameters with confidence interval.
#' @param object An object class of \code{poEMirtBoot}.
#' @param parameter A character or character vector, indicating what parameters ("alpha", "beta", "theta") you want to get. 
#' @param ci A float (< 1). Confidence interval.
#' @param ... Parameters to \code{summary()}.
#' @returns A tibble dataframe.
#' @importFrom dplyr %>% mutate select as_tibble bind_rows row_number bind_cols
#' @importFrom stringr str_remove str_c
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom stats qnorm
#' @export
summary.poEMirtBoot <- function(object, 
                                parameter = c('alpha', 'beta', 'theta'),
                                ci = 0.95,
                                ...) {
  cis <- c((1 - ci)/2, 1 - (1 - ci)/2)
  message('* Summarizing following parameters: ', paste(parameter, collapse = ', '), '.')
  out <- rep()
  if (any(parameter == 'alpha')) {
    if (!is.null(object$parameter$alpha)) {
      est_mean <- object$input$fit$parameter$alpha
      est_mean <- est_mean %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -.data$j,
          names_to = 'k',
          values_to = 'alpha'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          k = .data$k %>% 
            stringr::str_remove('y') %>% 
            as.integer(),
          index = stringr::str_c('[', .data$j, ', ', .data$k, ']'),
          reference = '[j, k]',
          parameter = 'alpha'
        ) %>% 
        dplyr::select(parameter, index, reference, estimate = alpha)
    
      sds <- object$standard_deviation$alpha
      sds <- sds %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -.data$j,
          names_to = 'k',
          values_to = 'sd'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(sd)
      
      out <- dplyr::bind_cols(est_mean, sds) %>% 
        dplyr::mutate(
          ci_lwr = .data$estimate + sd * stats::qnorm(cis[1]),
          ci_upr = .data$estimate + sd * stats::qnorm(cis[2])
        ) %>% 
        dplyr::bind_rows(out, .)
    } else {
      message('   - WARNING: You did not save `alpha` in your `poEMirtBoot` object. Set `control$save_item_parameters = TRUE` in `poEMirt_uncertainty` to obtain `alpha`.')
    }
  }
  if (any(parameter == 'beta')) {
    if (!is.null(object$parameter$beta)) {
      est_mean <- object$input$fit$parameter$beta
      est_mean <- est_mean %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -.data$j,
          names_to = 'k',
          values_to = 'beta'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          k = .data$k %>% 
            stringr::str_remove('y') %>% 
            as.integer(),
          index = stringr::str_c('[', .data$j, ', ', .data$k, ']'),
          reference = '[j, k]',
          parameter = 'beta'
        ) %>% 
        dplyr::select(parameter, index, reference, estimate = beta)
      
      sds <- object$standard_deviation$beta
      sds <- sds %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -.data$j,
          names_to = 'k',
          values_to = 'sd'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(.data$sd)
      
      out <- dplyr::bind_cols(est_mean, sds) %>% 
        dplyr::mutate(
          ci_lwr = .data$estimate + .data$sd * qnorm(cis[1]),
          ci_upr = .data$estimate + .data$sd * qnorm(cis[2])
        ) %>% 
        dplyr::bind_rows(out, .)
    } else {
      message('   - WARNING: You did not save `beta` in your `poEMirtBoot` object. Set `control$save_item_parameters = TRUE` in `poEMirt_uncertainty` to obtain `beta`.')
    }
  }
  if (any(parameter == 'theta')) {
    if (object$input$fit$info$model == 'static') {
      est_mean <- object$input$fit$parameter$theta
      est_mean <- est_mean %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(index = stringr::str_c('[', dplyr::row_number(), ']')) %>% 
        dplyr::mutate(
          reference = '[i]',
          parameter = 'theta'
        ) %>% 
        dplyr::select(parameter, index, reference, estimate = V1)
      
      sds <- object$standard_deviation$theta
      out <- dplyr::bind_cols(est_mean, sd = sds) %>%
        dplyr::mutate(
          ci_lwr = .data$estimate + .data$sd * qnorm(cis[1]),
          ci_upr = .data$estimate + .data$sd * qnorm(cis[2])
        ) %>% 
        dplyr::bind_rows(out, .)
    } else {
      est_mean <- object$input$fit$parameter$theta
      est_mean <- est_mean %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -.data$i,
          names_to = 't',
          values_to = 'theta'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          t = .data$t %>% 
            stringr::str_remove('V') %>% 
            as.integer(),
          index = stringr::str_c('[', .data$i, ', ', .data$t, ']'),
          reference = '[i, t]',
          parameter = 'theta'
        ) %>% 
        dplyr::select(.data$parameter, .data$index, .data$reference, estimate = .data$theta)
    
      sds <- object$standard_deviation$theta
      sds <- sds %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -.data$j,
          names_to = 'k',
          values_to = 'sd'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(.data$sd)
      out <- dplyr::bind_cols(est_mean, sd = sds) %>%
        dplyr::mutate(
          ci_lwr = .data$estimate + .data$sd * qnorm(cis[1]),
          ci_upr = .data$estimate + .data$sd * qnorm(cis[2])
        ) %>% 
        dplyr::bind_rows(out, .)
    }
  }
  return(out)
}

#' @title Summary function for \code{poEMirtGibbs}
#' @description \code{summary.poEMirtGibbs} returns a dataframe of estimated parameters with confidence interval.
#' @param object An object class of \code{poEMirtGibbs}.
#' @param parameter A character or character vector, indicating what parameters ("alpha", "beta", "theta") you want to get. 
#' @param ci A float (< 1). Confidence interval.
#' @param ... Parameters to \code{summary()}.
#' @returns A tibble dataframe.
#' @importFrom dplyr %>% mutate select as_tibble bind_rows row_number bind_cols
#' @importFrom stringr str_remove str_c
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom stats quantile median
#' @importFrom posterior rhat
#' @export
summary.poEMirtGibbs <- function(object, 
                                 parameter = c('alpha', 'beta', 'theta'),
                                 ci = 0.95,
                                 ...) {
  cis <- c((1 - ci)/2, 1 - (1 - ci)/2)
  message('* Summarizing following parameters: ', paste(parameter, collapse = ', '), '.')
  out <- rep()
  if (any(parameter == 'alpha')) {
    if (length(object$parameter$alpha) != 0) {
      est_mean <- dimstat(object$parameter$alpha, fun = 'mean', na.rm = TRUE)
      est_mean <- est_mean %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'alpha'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          k = k %>% 
            stringr::str_remove('y') %>% 
            as.integer(),
          index = stringr::str_c('[', j, ', ', k, ']'),
          reference = '[j, k]',
          parameter = 'alpha'
        ) %>% 
        dplyr::select(parameter, index, reference, mean = alpha)
      
      est_median <- apply(object$parameter$alpha, c(1, 2), stats::median, na.rm = TRUE)
      est_median <- est_median %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'median'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(median)
      
      sds <- object$standard_deviation$alpha
      sds <- sds %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'sd'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(sd)
      
      lwr <- apply(object$parameter$alpha, c(1, 2), quantile, probs = cis[1], na.rm = TRUE)
      upr <- apply(object$parameter$alpha, c(1, 2), quantile, probs = cis[2], na.rm = TRUE)
      lwr <- lwr %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'ci_lwr'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(ci_lwr)
      upr <- upr %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'ci_upr'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(ci_upr)
      
      rhats <- apply(object$parameter$alpha, c(1, 2), posterior::rhat)
      rhats <- rhats %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'rhat'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(rhat)
      out <- dplyr::bind_cols(est_mean, est_median, sds, lwr, upr, rhats) %>% 
        dplyr::bind_rows(out, .)
    } else {
      message('   - WARNING: You did not save `alpha` in your `poEMirtGibbs` object. Set `control$save_item_parameters = TRUE` in `poEMirt_uncertainty` to obtain `alpha`.')
    }
  }
  if (any(parameter == 'beta')) {
    if (length(object$parameter$beta) != 0) {
      est_mean <- dimstat(object$parameter$beta, fun = 'mean', na.rm = TRUE)
      est_mean <- est_mean %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'beta'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          k = k %>% 
            stringr::str_remove('y') %>% 
            as.integer(),
          index = stringr::str_c('[', j, ', ', k, ']'),
          reference = '[j, k]',
          parameter = 'beta'
        ) %>% 
        dplyr::select(parameter, index, reference, mean = beta)
      
      est_median <- apply(object$parameter$beta, c(1, 2), stats::median, na.rm = TRUE)
      est_median <- est_median %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'median'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(median)
      
      sds <- object$standard_deviation$beta
      sds <- sds %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'sd'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(sd)
      
      lwr <- apply(object$parameter$beta, c(1, 2), quantile, probs = cis[1], na.rm = TRUE)
      upr <- apply(object$parameter$beta, c(1, 2), quantile, probs = cis[2], na.rm = TRUE)
      lwr <- lwr %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'ci_lwr'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(ci_lwr)
      upr <- upr %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'ci_upr'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(ci_upr)
      
      rhats <- apply(object$parameter$beta, c(1, 2), posterior::rhat)
      rhats <- rhats %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'rhat'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(rhat)
      out <- dplyr::bind_cols(est_mean, est_median, sds, lwr, upr, rhats) %>% 
        dplyr::bind_rows(out, .)
    } else {
      message('   - WARNING: You did not save `beta` in your `poEMirtGibbs` object. Set `control$save_item_parameters = TRUE` in `poEMirt_uncertainty` to obtain `beta`.')
    }
  }
  if (any(parameter == 'theta')) {
    if (object$input$fit$info$model == 'static') {
      est_mean <- apply(object$parameter$theta, 1, mean, na.rm = TRUE)
      est_mean <- est_mean %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(index = stringr::str_c('[', dplyr::row_number(), ']')) %>% 
        dplyr::mutate(
          reference = '[i]',
          parameter = 'theta'
        ) %>% 
        dplyr::select(parameter, index, reference, mean = value)
      
      est_median <- apply(object$parameter$theta, 1, stats::median, na.rm = TRUE)
      sds <- object$standard_deviation$theta
      lwr <- apply(object$parameter$theta, 1, quantile, probs = cis[1], na.rm = TRUE)
      upr <- apply(object$parameter$theta, 1, quantile, probs = cis[2], na.rm = TRUE)
      rhats <- apply(object$parameter$theta, 1, posterior::rhat)
      out <- dplyr::bind_cols(est_mean, median = est_median, sd = sds, ci_lwr = lwr, ci_upr = upr, rhat = rhats) %>% 
        dplyr::bind_rows(out, .)
    } else {
      est_mean <- dimstat(object$parameter$theta, fun = 'mean', na.rm = TRUE)
      est_mean <- est_mean %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -i,
          names_to = 't',
          values_to = 'theta'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          t = t %>% 
            stringr::str_remove('V') %>% 
            as.integer(),
          index = stringr::str_c('[', i, ', ', t, ']'),
          reference = '[i, t]',
          parameter = 'theta'
        ) %>% 
        dplyr::select(parameter, index, reference, mean = theta)
      
      est_median <- apply(object$parameter$theta, c(1, 2), stats::median, na.rm = TRUE)
      est_median <- est_median %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'median'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(median)
      
      sds <- object$standard_deviation$theta
      sds[sds == 0] <- NA
      sds <- sds %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'sd'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(sd)
      
      lwr <- apply(object$parameter$theta, c(1, 2), quantile, probs = cis[1], na.rm = TRUE)
      upr <- apply(object$parameter$theta, c(1, 2), quantile, probs = cis[2], na.rm = TRUE)
      lwr <- lwr %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'ci_lwr'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(ci_lwr)
      upr <- upr %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'ci_upr'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(ci_upr)
      
      rhats <- apply(object$parameter$theta, c(1, 2), posterior::rhat)
      rhats <- rhats %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -j,
          names_to = 'k',
          values_to = 'rhat'
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(rhat)
      out <- dplyr::bind_cols(est_mean, est_median, sds, lwr, upr, rhats) %>% 
        dplyr::bind_rows(out, .)
    }
  }
  return(out)
}