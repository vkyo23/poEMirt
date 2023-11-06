#' @title Predict function
#' @description Predict function for `poEMirtFit`
#' @param object An object class of \code{poEMirtFit}
#' @param type A character, one of "prob" or "response". "prob" returns predicted probabilities of choices 
#' and "response" returns predicted value of Y.
#' @param seed An integer of random seed. This works when \code{type = "response"}.
#' @param ... Other arguments to \code{predict()}.
#' @return A dataframe of predicted values.
#' @rdname predict.poEMirtFit
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr %>% as_tibble mutate bind_cols select relocate
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_c str_replace
#' @importFrom rlang .data
#' @export

predict.poEMirtFit <- function(object, 
                               type = c("prob", "response"),
                               seed = NULL,
                               ...) 
{
  if (!is.null(seed)) set.seed(seed)
  type <- match.arg(type, c("prob", "response"))
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
  
  return(out)
}

