#' @title Predict function
#' @description Predict function for `poEMirtFit`
#' @param object An object class of \code{poEMirtFit}
#' @param type A character, one of "prob" or "response". "prob" returns predicted probabilities of choices 
#' and "response" returns predicted value of Y.
#' @param return Type of returning object. "array" or "dataframe".
#' @param seed An integer of random seed. This works when \code{type = "response"}.
#' @param ... Other arguments to \code{predict()}.
#' @return 3D array or dataframe
#' @rdname predict.poEMirtFit
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr %>% as_tibble mutate bind_cols select relocate
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_c str_replace
#' @importFrom rlang .data
#' @export

predict.poEMirtFit <- function(object, 
                               type = c("prob", "response"),
                               return = c("array", "dataframe"),
                               seed = NULL,
                               ...) {
  if (length(type) > 1) stop("`type` must be one of 'prob' or 'response'.")
  if (length(return) > 1) stop("`return` must be one of 'array' or 'dataframe'.")
  if (!is.null(seed)) set.seed(seed)
  type <- match.arg(type, c("prob", "response"))
  return <- match.arg(return, c("array", "dataframe"))
  object$parameter$alpha[is.na(object$parameter$alpha)] <- 0
  object$parameter$beta[is.na(object$parameter$beta)] <- 0
  object$parameter$theta[is.na(object$parameter$theta)] <- 0
  if (object$info$model != "static") {
    out <- prediction(
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
  if (type == "prob") {
    out[out == 0] <- NA
  }
  
  if (return == "dataframe") {
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

#' 
#' #' @title Predict function
#' #' @description Predict function for `poEMirtBoot
#' #' @param object An object class of \code{poEMirtBoot}
#' #' @param type A character, one of "prob" or "response". "prob" returns predicted probabilities of choices 
#' #' and "response" returns predicted value of Y.
#' #' @param return Type of returning object. "array" or "dataframe".
#' #' @param ... Other arguments to \code{predict()}.
#' #' @return A list of arrays or dataframe
#' #' @rdname predict.poEMirtBoot
#' #' @importFrom Rcpp sourceCpp
#' #' @importFrom dplyr %>% as_tibble 
#' #' @export
#' predict.poEMirtBoot <- function(object, 
#'                                 type = c("prob", "response"),
#'                                 return = c("array", "dataframe"),
#'                                 ...) {
#'   if (length(type) > 1) stop("`type` must be one of 'prob' or 'response'.")
#'   if (length(return) > 1) stop("`return` must be one of 'array' or 'dataframe'.")
#'   type <- match.arg(type, c("prob", "response"))
#'   return <- match.arg(return, c("array", "dataframe"))
#'   if (object$info$model != "static") {
#'     out <- prediction(
#'       N = object$info$data$data$trial,
#'       alpha = object$parameter$alpha,
#'       beta = object$parameter$beta,
#'       theta = object$parameter$theta,
#'       unique_categories = object$info$data$categories,
#'       item_timemap = object$info$data$dynamic$item_timemap,
#'       model = "dynamic",
#'       type = type
#'     )
#'   } else {
#'     out <- prediction(
#'       N = object$info$data$data$trial,
#'       alpha = object$parameter$alpha,
#'       beta = object$parameter$beta,
#'       theta = object$parameter$theta,
#'       unique_categories = object$info$data$categories,
#'       item_timemap = NA,
#'       model = "static",
#'       type = type
#'     )
#'   }
#'   dimnames(out) <- dimnames(object$info$data$data$response)
#'   return(out)
#' }
#' 
#' #' @title Predict function
#' #' @description Predict function for `poEMirtGibbs
#' #' @param object An object class of \code{poEMirtGibbs}
#' #' @param type A character, one of "prob" or "response". "prob" returns predicted probabilities of choices 
#' #' and "response" returns predicted value of Y.
#' #' @param ... Other arguments to \code{predict()}.
#' #' @return A dataframe
#' #' @rdname predict.poEMirtGibbs
#' #' @importFrom Rcpp sourceCpp
#' #' @export
#' predict.poEMirtGibbs <- function(object, 
#'                                  type = c("prob", "response"),
#'                                  ...) {
#'   if (length(type) > 1) stop("`type` must be one of 'prob' or 'response'.")
#'   type <- match.arg(type, c("prob", "response"))
#'   object$parameter$alpha[is.na(object$parameter$alpha)] <- 0
#'   object$parameter$beta[is.na(object$parameter$beta)] <- 0
#'   object$parameter$theta[is.na(object$parameter$theta)] <- 0
#'   B <- dim(object$parameter$alpha)[3]
#'   out <- array(NA, c(dim(object$input$fit$info$data$data$response), B))
#'   for (b in 1:B) {
#'     if (object$input$fit$info$model != "static") {
#'       out[, , , b] <- prediction(
#'         N = object$input$fit$info$data$data$trial,
#'         alpha = object$parameter$alpha[, , b],
#'         beta = object$parameter$beta[, , b],
#'         theta = object$parameter$theta[, , b],
#'         unique_categories = object$input$fit$info$data$categories,
#'         item_timemap = object$input$fit$info$data$dynamic$item_timemap,
#'         model = "dynamic",
#'         type = type
#'       )
#'     } else {
#'       out[, , , b] <- prediction(
#'         N = object$input$fit$info$data$data$trial,
#'         alpha = object$parameter$alpha[, , b],
#'         beta = object$parameter$beta[, , b],
#'         theta = object$parameter$theta[, b],
#'         unique_categories = object$input$fit$info$data$categories,
#'         item_timemap = NA,
#'         model = "static",
#'         type = type
#'       )
#'     }
#'   }
#'   dimnames(out) <- dimnames(object$input$fit$info$data$data$response)
#'   apply(out, c(1, 2, 3), mean)
#'   if (type == "prob") {
#'     out[out == 0] <- NA
#'   }
#'   
#'   dimnames(out) <- dimnames(object$info$data$data$response)
#'   return(out)
#' }