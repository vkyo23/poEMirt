#' @title Summary function for \code{poEMirtData}
#' @description \code{summary.poEMirtFit} returns a summary of \code{poEMirtData} object.
#' @param object An object class of \code{poEMirtFit}.
#' @param ... Other arguments to \code{summary()}.
#' @importFrom dplyr %>%
#' @export
summary.poEMirtData <- function(object, ...) {
  cat("----- poEMirtData summary -----\n")
  cat("* Data size:\n")
  cat("  - I = ", object$size$I, "\n")
  cat("  - J = ", object$size$J)
  if (exists("rep", object)) {
    al <- lapply(object$rep$processed, length) %>%
      unlist()
    cat(": Repeated j", paste0("[", sum(al > 1), " / ", length(al), "]"), "\n")
  } else {
    cat("\n")
  }
  if (exists("T", object$size)) {
    cat("  - T = ", object$size$T, "\n")
  }
  mi <- lapply(
    object$categories,
    max
  ) %>%
    unlist() %>%
    min()
  cat("  - min(Kj) = ", mi + 1, "\n")
  cat("  - max(Kj) = ", object$size$maxK, "\n")
  me <- lapply(
    object$categories,
    length
  ) %>% 
    unlist() %>%
    mean()
  cat("  - mean(Kj) =", round(me, 2), "\n")
  cat("* NA rate:", 
      paste0(round(sum(object$data$trial == 0) / length(object$data$trial), 3) * 100, "%"), "\n")
}

#' @title Summary function for \code{poEMirtFit}
#' @description \code{summary.poEMirtFit} returns a dataframe of estimated parameters.
#' @param object An object class of \code{poEMirtFit}.
#' @param parameter A character or character vector, indicating what parameters ("alpha", "beta", "theta") you want to get. 
#' @param ... Other arguments to \code{summary()}.
#' @returns A tibble dataframe.
#' @importFrom dplyr %>% mutate select as_tibble bind_rows row_number
#' @importFrom stringr str_c str_replace
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom rlang .data
#' @export
summary.poEMirtFit <- function(object, 
                               parameter = c("alpha", "beta", "theta"),
                               ...) {
  out <- rep()
  parameter <- match.arg(parameter, choices = c("alpha", "beta", "theta"), several.ok = TRUE)
  cat("* Summarizing following parameters:", paste(parameter, collapse = ", "), "\n")
  if (any(parameter == "alpha")) {
    tmp <- object$parameter$alpha %>% 
      dplyr::as_tibble(.name_repair = "unique") %>% 
      dplyr::mutate(j = rownames(object$parameter$alpha)) %>% 
      tidyr::pivot_longer(
        cols = -"j",
        names_to = "k",
        values_to = "alpha"
      ) %>% 
      tidyr::drop_na() %>% 
      dplyr::mutate(
        index = stringr::str_c("[", .data$j, ",", .data$k, "]"),
        reference = "[j,k]",
        parameter = "alpha"
      ) %>% 
      dplyr::select("parameter", "index", "reference", "estimate" = "alpha") 
    if (exists("rep", object$info$data)) {
      tmp <- tmp %>%
        dplyr::mutate(
          index = stringr::str_replace(.data$index, "-", ","),
          reference = "[t,j,k]"
        )
    }
    out <- dplyr::bind_rows(out, tmp)
  }
  if (any(parameter == "beta")) {
    tmp <- object$parameter$beta %>% 
      dplyr::as_tibble(.name_repair = "unique") %>% 
      dplyr::mutate(j = rownames(object$parameter$beta)) %>% 
      tidyr::pivot_longer(
        cols = -"j",
        names_to = "k",
        values_to = "beta"
      ) %>% 
      tidyr::drop_na() %>% 
      dplyr::mutate(
        index = stringr::str_c("[", .data$j, ",", .data$k, "]"),
        reference = "[j,k]",
        parameter = "beta"
      ) %>% 
      dplyr::select("parameter", "index", "reference", "estimate" = "beta") 
    if (exists("rep", object$info$data)) {
      tmp <- tmp %>%
        dplyr::mutate(
          index = stringr::str_replace(.data$index, "-", ","),
          reference = "[t,j,k]"
        )
    }
    out <- dplyr::bind_rows(out, tmp)
  }
  if (any(parameter == "theta")) {
    tmp <- object$parameter$theta %>% 
      dplyr::as_tibble(.name_repair = "unique")
    if (object$info$model == "static") {
      tmp <- tmp %>% 
        dplyr::mutate(index = stringr::str_c("[", rownames(object$parameter$theta), "]")) %>% 
        dplyr::mutate(
          reference = "[i]",
          parameter = "theta"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "estimate" = "1")
      out <- dplyr::bind_rows(out, tmp)
    } else {
      tmp <- tmp %>% 
        dplyr::mutate(i = rownames(object$parameter$theta)) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "theta"
        ) %>% 
        dplyr::mutate(
          t = .data$t %>% 
            as.integer(),
          index = stringr::str_c("[", .data$i, ",", .data$t, "]"),
          reference = "[i,t]",
          parameter = "theta"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "estimate" = "theta") %>% 
        filter(!is.na(.data$estimate))
      out <- dplyr::bind_rows(out, tmp)
    }
  }
  return(out)
}

#' @title Summary function for \code{poEMirtBoot}
#' @description \code{summary.poEMirtBoot} returns a dataframe of estimated parameters with confidence interval.
#' @param object An object class of \code{poEMirtBoot}.
#' @param parameter A character or character vector, indicating what parameters ("alpha", "beta", "theta") you want to get. 
#' @param ci A float (< 1). Confidence interval.
#' @param ... Other arguments to \code{summary()}.
#' @returns A tibble dataframe.
#' @importFrom dplyr %>% mutate select as_tibble bind_rows row_number bind_cols
#' @importFrom stringr str_c
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom stats qnorm
#' @importFrom rlang .data !! :=
#' @export
summary.poEMirtBoot <- function(object, 
                                parameter = c("alpha", "beta", "theta"),
                                ci = 0.95,
                                ...) {
  parameter <- match.arg(parameter, choices = c("alpha", "beta", "theta"), several.ok = TRUE)
  cis <- c((1 - ci)/2, 1 - (1 - ci)/2)
  cat("* Summarizing following parameters: ", paste(parameter, collapse = ", "), "\n")
  out <- rep()
  if (any(parameter == "alpha")) {
    if (!is.null(object$parameter$alpha)) {
      est_mean <- object$input$fit$parameter$alpha
      est_mean <- est_mean %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = rownames(object$parameter$alpha)) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "alpha"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          index = stringr::str_c("[", .data$j, ",", .data$k, "]"),
          reference = "[j,k]",
          parameter = "alpha"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "estimate" = "alpha")
      
      sds <- object$standard_deviation$alpha
      sds <- sds %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "sd"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("sd")
      
      tmp <- dplyr::bind_cols(est_mean, sds) %>% 
        dplyr::mutate(
          !!paste0("ci_lwr", ci * 100) := .data$estimate + sd * stats::qnorm(cis[1]),
          !!paste0("ci_upr", ci * 100) := .data$estimate + sd * stats::qnorm(cis[2])
        ) 
      if (exists("rep", object$input$fit$info$data)) {
        tmp <- tmp %>%
          dplyr::mutate(
            index = stringr::str_replace(.data$index, "-", ","),
            reference = "[t,j,k]"
          )
      }
      out <- dplyr::bind_rows(out, tmp)
    } else {
      warning("You did not save `alpha` in your `poEMirtBoot` object.\n Set `control$save_item_parameters = TRUE` in `poEMirt_uncertainty` to obtain `alpha`.")
    }
  }
  if (any(parameter == "beta")) {
    if (!is.null(object$parameter$beta)) {
      est_mean <- object$input$fit$parameter$beta
      est_mean <- est_mean %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = rownames(object$parameter$beta)) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "beta"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          index = stringr::str_c("[", .data$j, ",", .data$k, "]"),
          reference = "[j,k]",
          parameter = "beta"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "estimate" = "beta")
      
      sds <- object$standard_deviation$beta
      sds <- sds %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "sd"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("sd")
      
      tmp <- dplyr::bind_cols(est_mean, sds) %>% 
        dplyr::mutate(
          !!paste0("ci_lwr", ci * 100) := .data$estimate + sd * qnorm(cis[1]),
          !!paste0("ci_upr", ci * 100) := .data$estimate + sd * qnorm(cis[2])
        ) 
      if (exists("rep", object$input$fit$info$data)) {
        tmp <- tmp %>%
          dplyr::mutate(
            index = stringr::str_replace(.data$index, "-", ","),
            reference = "[t,j,k]"
          )
      }
      out <- dplyr::bind_rows(out, tmp)
    } else {
      warning("You did not save `beta` in your `poEMirtBoot` object.\n Set `control$save_item_parameters = TRUE` in `poEMirt_uncertainty` to obtain `beta`")
    }
  }
  if (any(parameter == "theta")) {
    if (object$input$fit$info$model == "static") {
      est_mean <- object$input$fit$parameter$theta
      est_mean <- est_mean %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(index = stringr::str_c("[", rownames(object$parameter$theta), "]")) %>% 
        dplyr::mutate(
          reference = "[i]",
          parameter = "theta"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "estimate" = "1")
      
      sds <- object$standard_deviation$theta
      tmp <- dplyr::bind_cols(est_mean, sd = sds) %>%
        dplyr::mutate(
          ci_lwr = .data$estimate + sd * qnorm(cis[1]),
          ci_upr = .data$estimate + sd * qnorm(cis[2])
        ) 
      out <- dplyr::bind_rows(out, tmp)
    } else {
      est_mean <- object$input$fit$parameter$theta
      est_mean <- est_mean %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(i = rownames(object$parameter$theta)) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "theta"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          t = .data$t %>% 
            as.integer(),
          index = stringr::str_c("[", .data$i, ",", .data$t, "]"),
          reference = "[i,t]",
          parameter = "theta"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "estimate" = "theta")
      
      sds <- object$standard_deviation$theta
      sds <- sds %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "sd"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("sd")
      tmp <- dplyr::bind_cols(est_mean, sd = sds) %>%
        dplyr::mutate(
          !!paste0("ci_lwr", ci * 100) := .data$estimate + sd * qnorm(cis[1]),
          !!paste0("ci_upr", ci * 100) := .data$estimate + sd * qnorm(cis[2])
        ) 
      out <- dplyr::bind_rows(out, tmp)
    }
  }
  return(out)
}

#' @title Summary function for \code{poEMirtGibbs}
#' @description \code{summary.poEMirtGibbs} returns a dataframe of estimated parameters with confidence interval.
#' @param object An object class of \code{poEMirtGibbs}.
#' @param parameter A character or character vector, indicating what parameters ("alpha", "beta", "theta") you want to get. 
#' @param ci A float (< 1). Credible interval.
#' @param ... Other arguments to \code{summary()}.
#' @returns A tibble dataframe.
#' @importFrom dplyr %>% mutate select as_tibble bind_rows row_number bind_cols
#' @importFrom stringr str_c str_replace
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom stats quantile median
#' @importFrom posterior rhat
#' @importFrom rlang .data
#' @export
summary.poEMirtGibbs <- function(object, 
                                 parameter = c("alpha", "beta", "theta"),
                                 ci = 0.95,
                                 ...) {
  parameter <- match.arg(parameter, choices = c("alpha", "beta", "theta"), several.ok = TRUE)
  cis <- c((1 - ci)/2, 1 - (1 - ci)/2)
  cat("* Summarizing following parameters:", paste(parameter, collapse = ", "), "\n")
  out <- rep()
  if (any(parameter == "alpha")) {
    if (length(object$parameter$alpha) != 0) {
      est_mean <- apply(object$parameter$alpha, c(1, 2), mean, na.rm = TRUE)
      est_mean <- est_mean %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = rownames(object$parameter$alpha)) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "alpha"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          index = stringr::str_c("[", .data$j, ",", .data$k, "]"),
          reference = "[j,k]",
          parameter = "alpha"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "mean" = "alpha")
      
      est_median <- apply(object$parameter$alpha, c(1, 2), stats::median, na.rm = TRUE)
      est_median <- est_median %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "median"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("median")
      
      sds <- object$standard_deviation$alpha
      sds <- sds %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "sd"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("sd")
      
      lwr <- apply(object$parameter$alpha, c(1, 2), stats::quantile, probs = cis[1], na.rm = TRUE)
      upr <- apply(object$parameter$alpha, c(1, 2), stats::quantile, probs = cis[2], na.rm = TRUE)
      lwr <- lwr %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "ci_lwr"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("ci_lwr")
      names(lwr) <- paste0(names(lwr), ci * 100)
      upr <- upr %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "ci_upr"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("ci_upr")
      names(upr) <- paste0(names(upr), ci * 100)
      
      rhats <- apply(object$parameter$alpha, c(1, 2), posterior::rhat)
      rhats <- rhats %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "rhat"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("rhat")
      tmp <- dplyr::bind_cols(est_mean, est_median, sds, lwr, upr, rhats)  
      if (exists("rep", object$input$fit$info$data)) {
        tmp <- tmp %>%
          dplyr::mutate(
            index = stringr::str_replace(.data$index, "-", ","),
            reference = "[t,j,k]"
          )
      }
      out <- dplyr::bind_rows(out, tmp)
    } else {
      warning("You did not save `alpha` in your `poEMirtGibbs` object.\n Set `control$save_item_parameters = TRUE` in `poEMirt_uncertainty` to obtain `alpha`.")
    }
  }
  if (any(parameter == "beta")) {
    if (length(object$parameter$beta) != 0) {
      est_mean <- apply(object$parameter$beta, c(1, 2), mean, na.rm = TRUE)
      est_mean <- est_mean %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = rownames(object$parameter$beta)) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "beta"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          index = stringr::str_c("[", .data$j, ",", .data$k, "]"),
          reference = "[j,k]",
          parameter = "beta"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "mean" = "beta")
      
      est_median <- apply(object$parameter$beta, c(1, 2), stats::median, na.rm = TRUE)
      est_median <- est_median %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "median"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("median")
      
      sds <- object$standard_deviation$beta
      sds <- sds %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "sd"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("sd")
      
      lwr <- apply(object$parameter$beta, c(1, 2), stats::quantile, probs = cis[1], na.rm = TRUE)
      upr <- apply(object$parameter$beta, c(1, 2), stats::quantile, probs = cis[2], na.rm = TRUE)
      lwr <- lwr %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "ci_lwr"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("ci_lwr")
      names(lwr) <- paste0(names(lwr), ci * 100)
      upr <- upr %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "ci_upr"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("ci_upr")
      names(upr) <- paste0(names(upr), ci * 100)
      
      rhats <- apply(object$parameter$beta, c(1, 2), posterior::rhat)
      rhats <- rhats %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(j = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"j",
          names_to = "k",
          values_to = "rhat"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("rhat")
      tmp <- dplyr::bind_cols(est_mean, est_median, sds, lwr, upr, rhats) 
      if (exists("rep", object$input$fit$info$data)) {
        tmp <- tmp %>%
          dplyr::mutate(
            index = stringr::str_replace(.data$index, "-", ","),
            reference = "[t,j,k]"
          )
      }
      out <- dplyr::bind_rows(out, tmp)
    } else {
      warning("You did not save `beta` in your `poEMirtGibbs` object.\n Set `control$save_item_parameters = TRUE` in `poEMirt_uncertainty` to obtain `beta`.")
    }
  }
  if (any(parameter == "theta")) {
    if (object$input$fit$info$model == "static") {
      est_mean <- apply(object$parameter$theta, 1, mean, na.rm = TRUE)
      est_mean <- est_mean %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(index = stringr::str_c("[", rownames(object$parameter$theta), "]")) %>% 
        dplyr::mutate(
          reference = "[i]",
          parameter = "theta"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "mean" = "value")
      
      est_median <- apply(object$parameter$theta, 1, stats::median, na.rm = TRUE)
      sds <- object$standard_deviation$theta
      lwr <- apply(object$parameter$theta, 1, stats::quantile, probs = cis[1], na.rm = TRUE)
      upr <- apply(object$parameter$theta, 1, stats::quantile, probs = cis[2], na.rm = TRUE)
      rhats <- apply(object$parameter$theta, 1, posterior::rhat)
      tmp <- dplyr::bind_cols(est_mean, median = est_median, sd = sds, ci_lwr = lwr, ci_upr = upr, rhat = rhats) 
      out <- dplyr::bind_rows(out, tmp)
    } else {
      est_mean <- apply(object$parameter$theta, c(1, 2), mean, na.rm = TRUE)
      est_mean <- est_mean %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(i = rownames(object$parameter$theta)) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "theta"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(
          t = .data$t %>% 
            as.integer(),
          index = stringr::str_c("[", .data$i, ",", .data$t, "]"),
          reference = "[i,t]",
          parameter = "theta"
        ) %>% 
        dplyr::select("parameter", "index", "reference", "mean" = "theta")
      
      est_median <- apply(object$parameter$theta, c(1, 2), stats::median, na.rm = TRUE)
      est_median <- est_median %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "median"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("median")
      
      sds <- object$standard_deviation$theta
      sds[sds == 0] <- NA
      sds <- sds %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "sd"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("sd")
      
      lwr <- apply(object$parameter$theta, c(1, 2), stats::quantile, probs = cis[1], na.rm = TRUE)
      upr <- apply(object$parameter$theta, c(1, 2), stats::quantile, probs = cis[2], na.rm = TRUE)
      lwr <- lwr %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "ci_lwr"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("ci_lwr")
      names(lwr) <- paste0(names(lwr), ci * 100)
      upr <- upr %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "ci_upr"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("ci_upr")
      names(upr) <- paste0(names(upr), ci * 100)
      rhats <- apply(object$parameter$theta, c(1, 2), posterior::rhat)
      rhats <- rhats %>% 
        dplyr::as_tibble(.name_repair = "unique") %>% 
        dplyr::mutate(i = dplyr::row_number()) %>% 
        tidyr::pivot_longer(
          cols = -"i",
          names_to = "t",
          values_to = "rhat"
        ) %>% 
        tidyr::drop_na() %>% 
        dplyr::select("rhat")
      tmp <- dplyr::bind_cols(est_mean, est_median, sds, lwr, upr, rhats) 
      out <- dplyr::bind_rows(out, tmp)
    }
  }
  return(out)
}