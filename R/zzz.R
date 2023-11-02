.onAttach <- function(...) {
  packageStartupMessage("* `poEMirt`: A fast item response theory models for public opinion data analysis")
}

.onUnload <- function(libpath) {
  library.dynam.unload("poEMirt", libpath)
}