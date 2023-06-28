
#' @title Manages a Sequence of Seeds
#'
#' @description Returns a series of random number seeds from a flexible set of
#' inputs.
#'
#' @param ref.seed a value that will be coerced to integer OR an object with
#' length equal to `n` OR `NULL`: either, respectively
#'  - the starting point for the sequence of seeds
#'  - the sequence of seeds
#'  - the indicator to *not* set seeds
#'
#' @param n a positive integer: the number of seeds to return
#'
#' @return depending on the input format:
#'  - a sequence of seeds, `ref.seed + 1L:n - 1L`
#'  - the same sequence of seeds provided
#'  - a list of length `n`, with all `NULL` elements
#'
#' @examples
#' seq_seeds # to view function internals
#' tenseeds <- seq_seeds(0, 10)
#' tenseeds
#' tenmoreseeds <- seq_seeds(0:9, 10)
#' tenmoreseeds
#' # note, won't work with nonsense input:
#' try(seq_seed(0:10, 10))
#' # intended for use in `apply`-family invocations of functions that optionally
#' # want a seed:
#' seq_seed(n=10)
#'
#' @export
seq_seeds <- function(ref.seed, n) {
  n |> check_integer() |> check_positive() |> check_scalar()
  if (!missing(ref.seed)) {
    if (length(ref.seed) != 1) {
      length(ref.seed) |> check_equal(n)
      seeds <- ref.seed
    } else {
      seeds <- ref.seed + seq_len(n) - 1L
    }
  } else {
    seeds <- rep(list(NULL), times = n)
  }
  return(seeds)
}

#' @title Manages a Sequence of Parameters
#'
#' @description Returns a series of parameters from a flexible set of
#' inputs.
#'
#' @param ref.seed a named list OR an unnamed list of length equal to `n`:
#' respectively
#'  - the parameter values to be replicated into a sequence
#'  - the sequence of parameter values
#'
#' @param n a positive integer: the number of parameter combinations to return
#'
#' @return depending on the input format:
#'  - a list of clones of parameters
#'  - the same sequence of parameters provided
#'
#' @export
#' @family utilities
seq_parms <- function(ref.parms, n) {
  n |> check_integer() |> check_positive() |> check_scalar()
  if (is.null(names(ref.parms))) {
    length(ref.parms) |> check_equal(n)
    parms <- ref.parms
  } else {
    parms <- rep(list(ref.parms), n)
  }
  return(parms)
}

#' @title Sequence from Range
#'
#' @description
#' Returns as sequence between two ends as returned by [range()].
#'
#' @param r an output from [range()]
#'
#' @inheritDotParams base::range
#'
#' @export
#' @family utilities
seq_between <- function(r, ...) seq(from = r[1], to = r[2], ...)

#' @title Get Year
#'
#' @description
#' Gets the year from a date (the current date, by default).
#'
#' @param date a object coerce-able via [as.Date()]
#'
#' @return an integer, the 4 digit year corresponding to `date`
#'
#' @export
#' @family utilities
theyear <- function(date = Sys.Date()) {
  date |> as.Date() |> trunc("years") |> format("%Y") |> as.integer()
}

#' @title Launch an ICI3D Tutorial
#'
#' @param tutorial a string: the tutorial name (must be among the defaults)
#'
#' @export
launch <- function(tutorial = c("dynfever")) {
  tutorial <- match.arg(tutorial)
  learnr::run_tutorial(name = tutorial, package = "ICI3D")
}
