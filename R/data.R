#' @title Simulated data (dynamic)
#' @format A tbl_df with 200 individuals, 400 items, 40 time-periods (10 items per time-period), and max(K_j) = 5. 
#'  \describe{
#'   \item{i}{Individual ID}
#'   \item{j}{Item ID}
#'   \item{t}{Time-period}
#'   \item{y1}{Response to k = 1}
#'   \item{y2}{Response to k = 2}
#'   \item{y3}{Response to k = 3}
#'   \item{y4}{Response to k = 4}
#'   \item{y5}{Response to k = 5}
#' }
"sim_data_dynamic"

#' @title Simulated data (static)
#' @format A tbl_df with 500 individuals, 100 items, and max(K_j) = 5. 
#'  \describe{
#'   \item{i}{Individual ID}
#'   \item{j}{Item ID}
#'   \item{y1}{Response to k = 1}
#'   \item{y2}{Response to k = 2}
#'   \item{y3}{Response to k = 3}
#'   \item{y4}{Response to k = 4}
#'   \item{y5}{Response to k = 5}
#' }
"sim_data_static"