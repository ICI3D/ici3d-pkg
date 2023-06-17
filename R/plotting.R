
# geom_grid <- function(
#     ys, FUN = as.character
# ) {
#   geom_texthline(
#     aes(yintercept = ys, label = label),
#     data = data.table(ys = ys, label = FUN(ys)), inherit.aes = FALSE,
#     gap = TRUE, hjust = 0.325, color = "grey85", fontface = "bold"
#   )
# }

#' @title Extract Positions to Place Labels on Grid Lines
#'
#' @description Intended for use with [ggplot2] `geom`s `data` argument, uses
#' `data.table`-like syntax to identify label positions + values.
#'
#' @param dt a `data.table`: generally intended to be the data for the whole
#' `ggplot` object.
#'
#' @param axiscol a bare variable name: the column in the `data.table` which the
#' scale will be over.
#'
#' @param scalefun a function: should accept a vector of the same type as
#' represented by `axiscol` and return a vector of that type.
#'
#' @param labelfun a function: should accept a vector of the same type as
#' represented by `axiscol` and return a `character` vector (of labels).
#'
#' @param by optional expression of the same type as [data.table::data.table]
#' `by` argument.
#'
#' @param ... other named arguments: these will be added as columns to the
#' returned `data.table`
#'
#' @import scales data.table
#' @export
data_grid_labels <- function(
  dt, axiscol,
  scalefun = scales::pretty_breaks(), labelfun = as.character,
  by, ..., exclude.zero = TRUE
) {
  res <- eval(substitute(dt[,{ # NSE to let us use bare names in arguments
    tmp <- scalefun(axiscol)
    if (length(tmp) > 2) { # might want to exclude min/max labels
      mm <- range(axiscol)
      if (mm[1] > (tmp[1]-tmp[2])/2 + tmp[2]) {
        tmp <- tail(tmp, -1)
      }
      if (mm[2] < (tail(tmp,1) - diff(tail(tmp, 2))/2)) {
        tmp <- head(tmp, -1)
      }
    }
    if (exclude.zero) tmp <- tmp[tmp != 0]
    .(tmp = tmp)
  }, by=by]))
  res |> setnames("tmp", substitute(axiscol) |> as.character())
  adds <- list(...)
  if (length(adds)) {
    res[, names(adds) := adds ]
  }
  sc <- substitute(axiscol)
  res[, label := eval(sc) |> labelfun() ]
  return(res)
}
