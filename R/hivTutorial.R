#' @import deSolve
NULL

#' @export
hivTutorial <- function() rmarkdown::run(system.file("hiv", "index.Rmd", package = "ICI3D"))
