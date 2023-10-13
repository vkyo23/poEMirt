#' Predict function for `poEMirtFit`
#' @param object An object class of \code{poEMirtFit}
#' @param ... Auxiliary parameters to \code{predict()}.
#' @param type A character, one of "prob" or "response". "prob" returns predicted probabilities of choices 
#' and "response" returns predicted value of Y.
#' @rdname predict.poEMirtFit
#' @importFrom Rcpp sourceCpp
#' @export

predict.poEMirtFit <- function(object, 
                               type = c('prob', 'response'),
                               ...) {
  type <- match.arg(type, c('prob', 'response'))
  if (object$info$model != 'static') {
    out <- prediction(
      N = object$info$data$trial,
      alpha = object$parameter$alpha,
      beta = object$parameter$beta,
      theta = object$parameter$theta,
      unique_categories = object$info$data$categories,
      item_timemap = object$info$data$dynamic$item_timemap - 1,
      model = 'dynamic',
      type = type
    )
  } else {
    out <- prediction(
      N = object$info$data$trial,
      alpha = object$parameter$alpha,
      beta = object$parameter$beta,
      theta = object$parameter$theta,
      unique_categories = object$info$data$categories,
      item_timemap = NA,
      model = 'static',
      type = type
    )
  }
  dimnames(out) <- dimnames(object$info$data$response)
  return(out)
}
