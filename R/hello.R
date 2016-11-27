# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}

#' @export
dynamicalFever <- function() runApp(system.file('dynFev', package='DAIDD'))

#' @export
reedFrost <- function() runApp(system.file('reedFrost', package='DAIDD'))

#' @export
heterogenietyTutorial <- function() runApp(system.file('hetTut', package='DAIDD'))

#' @import shiny
NULL

#' @import ggplot2
NULL

#' @import data.table
NULL
