app.ui.generic <- function(title, ...) shinyUI(fluidPage(
  titlePanel(title), navlistPanel(..., widths = c(2, 10))
))
