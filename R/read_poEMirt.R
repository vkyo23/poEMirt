#' @title Preparing for \code{poEMirt()}
#' @description This function creates an input data for \code{poEMirt} models.
#' 
#' @param dataframe A data.frame or tbl_df object.
#' @param responses Names of columns of responses. Last response category for each j will be the baseline.
#' @param i A name of individual IDs column. Must be an integer column.
#' @param j A name of item IDs column. Must be an integer column.
#' @param t A name of time period IDs column (for `dynamic` model). Must be an integer column.
#' @param smooth A character, the smoothing method. 
#' \itemize{
#'   \item \code{"no"} No smoothing.
#'   \item \code{"default"} Using from the start time of i to T.
#'   \item \code{"finite"} Using from the start time of it to the end time.
#'  }
#' 
#' @return A \code{poEMirtData} object.
#' 
#' @importFrom rlang sym !! enquo .data :=
#' @importFrom dplyr %>% distinct select filter right_join left_join arrange mutate across everything all_of group_by summarize summarize_at mutate_at row_number
#' @importFrom tidyr expand_grid spread 
#' @importFrom stringr str_c str_split
#' @importFrom Rcpp sourceCpp
#' @export
#' 
#' @examples 
#' \dontrun{
#' data("sim_data_dynamic")
#' 
#' # Convert into poEMirt-readable data
#' data <- read_poEMirt(
#'   dataframe = sim_data_dynamic,
#'   responses = paste0('y', 1:5),
#'   i = "i",
#'   j = "j",
#'   t = "t"
#' )
#' }

read_poEMirt <- function(dataframe, 
                         responses, 
                         i, 
                         j, 
                         t = NULL, 
                         smooth = NULL) {
  # Input check
  if (!class(dataframe)[1] %in% c("tbl_df", "data.frame")) stop("`dataframe` should be a data.frame or tbl_df.")
  if (is.null(smooth)) smooth <- "default"
  smooth <- match.arg(smooth, choices = c("no", "default", "finite"))
  if (length(smooth) > 1) stop("`smooth` must be one of 'no', 'default', or 'finite'.")
  
  # Size
  i <- rlang::sym(i)
  j <- rlang::sym(j)
  responses <- rlang::enquo(responses)
  if (!is.null(t)) {
    t <- rlang::sym(t)
    # Check whether j is repeated id or not
    tmp <- dataframe %>% 
      dplyr::distinct(!!j, !!t) %>% 
      dplyr::group_by(!!j) %>% 
      dplyr::summarize(count = dplyr::n())
    if (sum(tmp$count) > nrow(tmp)) {
      repeated_j <- TRUE
      cat("* Detecting repeated j", paste0("[", nrow(tmp[tmp$count >= 2, ]), " / ", nrow(tmp), "]"), ":\n")
      cat("  - Set `alpha_fix = TRUE` in `poEMirt()` to fix alpha of same repeated items.\n")
      
      # Create unique item ids
      j_val <- tmp %>% 
        nrow()
      j_val <- 10^nchar(j_val)
      t_val <- dataframe %>% 
        dplyr::distinct(!!t) %>% 
        nrow()
      t_val <- 10^nchar(t_val)
      dataframe <- dataframe %>% 
        dplyr::mutate(
          unique_j = stringr::str_c(t_val + !!t, "-", j_val + !!j)
        )
      
      # Identifying repeated j
      rep_j <- dataframe %>% 
        dplyr::distinct(.data$unique_j, !!j) %>% 
        dplyr::arrange(.data$unique_j) 
      
      # Replace j with unqiue_j
      dataframe <- dataframe %>% 
        dplyr::mutate(!!j := .data$unique_j)
    } else {
      repeated_j <- FALSE
    }
  }
  
  # Exclude some items
  excl <- dataframe %>% 
    dplyr::group_by(!!j) %>% 
    dplyr::summarize_at(
      .vars = dplyr::vars(dplyr::all_of(!!responses)),
      .funs = function(x) sum(x, na.rm = TRUE)
    ) 
  jname <- as.vector(as.matrix(excl[, 1]))
  excl <- excl %>% 
    dplyr::select(-!!j)
  excl <- excl %>% 
    dplyr::mutate(
      all = rowSums(excl)
    ) 
  excl <- excl %>% 
    dplyr::mutate_at(
      .vars = dplyr::vars(dplyr::all_of(!!responses)),
      .funs = function(x) x == excl$all
    ) %>% 
    dplyr::select(-"all") 
  excl <- excl %>% 
    dplyr::mutate(
      exc = rowSums(excl)
    ) %>% 
    dplyr::mutate(
      j = jname
    ) %>% 
    dplyr::filter(.data$exc == 1)
  if (nrow(excl) > 0) {
    cat("* Remove following items due to no variation in responses\n")
    cat("  -", excl$j, "\n")
    dataframe <- dataframe %>% 
      dplyr::filter(
        !eval(j) %in% excl$j
      )
    if (!is.null(t)) {
      if (repeated_j) {
        rep_j <- rep_j %>% 
          dplyr::filter(
            !eval(j) %in% excl$j
          )
      }
    }
  }
  
  if (!is.null(t)) {
    # Dynamic 
    dataframe <- dataframe %>% 
      dplyr::select(i = !!i, j = !!j, t = !!t, dplyr::all_of(!!responses))
    T <- dataframe %>% 
      dplyr::distinct(t) %>% 
      nrow()
    item_timemap <- dataframe %>% 
      dplyr::distinct(j, t) %>% 
      dplyr::arrange(j)
  } else {
    dataframe <- dataframe %>% 
      dplyr::select(i = !!i, j = !!j, dplyr::all_of(!!responses))
  }
  I <- dataframe %>% 
    dplyr::distinct(i) %>% 
    nrow()
  J <- dataframe %>% 
    dplyr::distinct(j) %>% 
    nrow()
  maxK <- dataframe %>% 
    dplyr::select(dplyr::all_of(!!responses)) %>% 
    ncol()
  
  # Response array
  Y <- array(NA, dim = c(I, J, maxK))
  Ks <- dataframe %>% 
    dplyr::select(dplyr::all_of(!!responses)) %>% 
    colnames()
  for (kk in 1:length(Ks)) {
    tmp <- dataframe %>% 
      dplyr::select(i, j, y = dplyr::all_of(Ks[kk])) %>% 
      tidyr::spread(key = j, value = .data$y) %>%
      arrange(i)
    Y[, , kk] <- tmp[, -1] %>%
      as.matrix() 
  }
  
  rownames(Y) <- sort(unique(dataframe$i))
  coln <- sort(unique(dataframe$j))
  if (!is.null(t)) {
    if (repeated_j) {
      tmp1 <- stringr::str_split(coln, "-")
      tt <- lapply(
        tmp1,
        function(x) {
          as.integer(x[1]) - t_val
        }
      ) %>% 
        unlist()
      jj <- lapply(
        tmp1,
        function(x) {
          as.integer(x[2]) - j_val
        }
      ) %>% 
        unlist()
      coln <- paste0(tt, "-", jj)
    }
  }
  colnames(Y) <- coln
  dimnames(Y)[[3]] <- Ks
  
  # Categories
  categories <- apply(colSums(Y, na.rm = TRUE), 1, function(x) which(x != 0)-1, simplify = FALSE)
  
  # Number of trials
  n <- apply(Y, c(1, 2), sum, na.rm = TRUE)
  
  # Get data information
  sb <- construct_sb_auxs(Y, n, categories)
  
  # Create an output list
  L <- list(
    data = list(
      response = Y,
      trial = n, 
      modelinput = sb
    ),
    size = list(
      I = I,
      J = J,
      maxK = maxK
    ),
    categories = categories
  )
  
  if (!is.null(t)) {
    # Repeated j
    if (repeated_j) {
      rep_j_mat <- rep_j %>% 
        dplyr::mutate(dum = 1) %>% 
        tidyr::spread(
          key = .data$unique_j,
          value = .data$dum
        )
      rep_j_mat <- as.matrix(rep_j_mat[, -1])
      rep_j_mat[is.na(rep_j_mat)] <- 0
      prc <- apply(rep_j_mat, 1, function(x) which(x == 1) - 1, simplify = FALSE)
      
      # Get index
      L$rep <- list(
        processed = prc,
        raw = rep_j_mat,
        flag = repeated_j
      )
    }
    
    # Time-map
    tmp <- dataframe %>% 
      dplyr::distinct(i, t) %>% 
      dplyr::arrange(i) %>% 
      dplyr::mutate(dum = 1) %>% 
      tidyr::spread(key = t, value = .data$dum) 
    timemap <- tmp[, -1] %>%
      dplyr::mutate(
        dplyr::across(
          .cols = dplyr::everything(),
          .fns = function(x) ifelse(is.na(x), 0, x)
        )
      ) %>% 
      as.matrix()
    rownames(timemap) <- rownames(Y)
    timemap2 <- timemap
    if (smooth == "default") {
      for (ii in 1:I) {
        Ti_start <- names(which.min(which(timemap[ii, ] == 1)))
        timemap[ii, Ti_start:T] <- 1
      }
    } else if (smooth == "finite") {
      for (ii in 1:I) {
        Ti_start <- names(which.min(which(timemap[ii, ] == 1)))
        Ti_end <- names(which.max(which(timemap[ii, ] == 1)))
        timemap[ii, Ti_start:Ti_end] <- 1
      }
    }
    
    # Item time-map
    item_timemap <- item_timemap$t - 1
    names(item_timemap) <- colnames(Y)
    
    L$size$T <- T
    
    index <- get_dynamic_info(n, item_timemap, timemap)
    
    L$dynamic <- list(
      timemap = timemap,
      timemap2 = timemap2,
      item_timemap = item_timemap,
      smooth = smooth,
      index = index
    )
  }
  
  class(L) <- c("poEMirtData", class(L))
  
  return(L)
}
