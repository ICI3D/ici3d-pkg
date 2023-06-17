
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
      seeds <- ref.seed + (1L:n) - 1L
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
#'  - the paramaeter values to be replicated into a sequence
#'  - the sequence of parameter values
#'
#' @param n a positive integer: the number of parameter combinations to return
#'
#' @return depending on the input format:
#'  - a list of clones parameters
#'  - the same sequence of parameters provided
#'
#' @export
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

#' @export
seq_between <- function(r, ...) seq(from = r[1], to = r[2], ...)

#' @export
theyear <- function(date = Sys.Date()) {
  date |> as.Date() |> trunc("years") |> format("%Y") |> as.integer()
}
