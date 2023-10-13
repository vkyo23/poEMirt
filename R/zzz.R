.onAttach <- function(...) {
  d <- date()
  tmp <- regexpr("[0-9]{4}", d)
  year <- substr(d, tmp[1], tmp[1] + attr(tmp, "match.length") - 1)
  packageStartupMessage("* poEMirt (", year, "): A fast item response theory models for public opinion data analysis")
}
