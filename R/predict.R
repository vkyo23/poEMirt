#' @title Predict function
#' @description Predict function for `poEMirtFit`
#' @param object An object class of \code{poEMirtFit}
#' @return A dataframe of predicted values.
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr %>% as_tibble mutate bind_cols select relocate
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_c str_replace
#' @importFrom rlang .data
#' @useDynLib poEMirt, .registration = TRUE
#' @noRd
#' @export

predict.poEMirtFit <- function(object, 
                               ...) 
{
  predict1(object, ...)
}
