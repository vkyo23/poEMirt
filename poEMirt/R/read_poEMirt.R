#' @title Preparing for \code{poEMirt()}
#' @description This function creates an input data for \code{poEMirt} models.
#' 
#' @param data A data.frame or tbl_df object.
#' @param responses Names of columns of responses. Last response category will be the baseline.
#' @param i A name of individual IDs column.
#' @param j A name of item IDs column.
#' @param t A name of time period IDs column (for `dynamic` model).
#' @param no_smooth A boot. If TRUE, the model does not estimate latent traits of i in missing years (for `dynamic` model). Default is FALSE.
#' @param matched_j A name of column that contains item IDs of matched item.
#' @return A \code{poEMirtData} object.
#' 
#' @importFrom rlang sym !! enquo
#' @importFrom dplyr %>% distinct select filter right_join left_join arrange mutate across everything all_of group_by summarise_at mutate_at row_number
#' @importFrom tidyr expand_grid spread 
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
#' }

read_poEMirt <- function(data, 
                         responses, 
                         i, 
                         j, 
                         t = NULL, 
                         no_smooth = FALSE, 
                         matched_j = NULL) {
  # Input check
  if (!class(data)[1] %in% c('tbl_df', 'data.frame')) stop('`data` should be a data.frame or tbl_df.')
  
  # Size
  i <- rlang::sym(i)
  j <- rlang::sym(j)
  responses <- rlang::enquo(responses)
  
  # Exclude some items
  excl <- data %>% 
    dplyr::group_by(!!j) %>% 
    dplyr::summarise_at(
      .vars = dplyr::vars(dplyr::all_of(!!responses)),
      .funs = function(x) sum(x, na.rm = TRUE)
    ) 
  jname <- excl$item_year_id
  excl <- excl |> 
    dplyr::select(-!!j) %>% 
    dplyr::mutate(
      all = rowSums(.)
    ) 
  excl <- excl %>% 
    dplyr::mutate_at(
      .vars = dplyr::vars(dplyr::all_of(!!responses)),
      .funs = function(x) x == excl$all
    ) %>% 
    dplyr::select(-all) %>% 
    dplyr::mutate(
      exc = rowSums(.)
    ) %>% 
    dplyr::mutate(
      j = jname
    ) %>% 
    dplyr::filter(exc == 1)
  if (nrow(excl) > 0) {
    cat('* Remove following items due to no variation in responses\n')
    cat('  -', excl$j, '\n')
    data <- data %>% 
      dplyr::filter(
        !eval(j) %in% excl$j
      )
  }
  
  if (!is.null(t)) {
    # Dynamic 
    t <- rlang::sym(t)
    
    if (!is.null(matched_j)) {
      matched_j <- rlang::sym(matched_j)
      data <- data %>% 
        dplyr::select(i = !!i, j = !!j, mj = !!matched_j, t = !!t,dplyr::all_of(!!responses))
    } else {
      data <- data %>% 
        dplyr::select(i = !!i, j = !!j, t = !!t, dplyr::all_of(!!responses))
    }
    
    T <- data %>% 
      dplyr::distinct(t) %>% 
      nrow()
    item_timemap <- data %>% 
      dplyr::distinct(j, t) %>% 
      dplyr::arrange(j)
  } else {
    data <- data %>% 
      dplyr::select(i = !!i, j = !!j, dplyr::all_of(!!responses))
  }
  I <- data %>% 
    dplyr::distinct(i) %>% 
    nrow()
  J <- data %>% 
    dplyr::distinct(j) %>% 
    nrow()
  maxK <- data %>% 
    dplyr::select(dplyr::all_of(!!responses)) %>% 
    ncol()
  
  # Response array
  Y <- array(NA, dim = c(I, J, maxK))
  Ks <- data %>% 
    dplyr::select(dplyr::all_of(!!responses)) %>% 
    colnames()
  for (kk in 1:length(Ks)) {
    Y[, , kk] <- data %>% 
      dplyr::select(i, j, y = dplyr::all_of(Ks[kk])) %>% 
      tidyr::spread(key = j, value = y) %>%
      arrange(i) %>%
      .[, -1] %>%
      as.matrix() 
  }
  
  rownames(Y) <- sort(unique(data$i))
  colnames(Y) <- sort(unique(data$j))
  dimnames(Y)[[3]] <- Ks
  
  # Categories
  categories <- apply(colSums(Y, na.rm = TRUE), 1, function(x) which(x != 0), simplify = FALSE)
  
  # Number of trials
  n <- apply(Y, c(1, 2), sum, na.rm = TRUE)
  
  # Create an output list
  L <- list(
    response = Y, 
    trial = n, 
    size = list(
      I = I,
      J = J,
      maxK = maxK
    ),
    categories = categories
  )
  
  if (!is.null(t)) {
    # Time-map
    timemap <- data %>% 
      dplyr::distinct(i, t) %>% 
      dplyr::arrange(i) %>% 
      dplyr::mutate(dum = 1) %>% 
      tidyr::spread(key = t, value = dum) %>% 
      .[, -1] %>% 
      dplyr::mutate(
        dplyr::across(
          .cols = dplyr::everything(),
          .fns = function(x) ifelse(is.na(x), 0, x)
        )
      ) %>% 
      as.matrix()
    rownames(timemap) <- rownames(Y)
    timemap2 <- timemap
    if (!no_smooth) {
      for (ii in 1:I) {
        Ti_start <- names(which.min(which(timemap[ii, ] == 1)))
        timemap[ii, Ti_start:T] <- 1
      }
    }
    
    # Item time-map
    item_timemap <- item_timemap$t
    names(item_timemap) <- colnames(Y)
    
    # Matched items
    if (!is.null(matched_j)) {
      item_match <- data %>% 
        dplyr::distinct(j, mj) %>% 
        dplyr::arrange(j) %>% 
        .$j
    } else {
      item_match <- rep(NA, J)
    }
    
    L$size$T <- T
    L$dynamic <- list(
      timemap = timemap,
      timemap2 = timemap2,
      item_timemap = item_timemap,
      item_match = item_match,
      no_smooth = no_smooth
    )
  }
  
  class(L) <- c('poEMirtData', class(L))
  return(L)
}
