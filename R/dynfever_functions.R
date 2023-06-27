
#' @title Initialize a Dynamical Fever Population
#'
#' @param dogs an integer: how many dogs in the population?
#'
#' @param humans an integer: how many humans in the population?
#'
#' @return a named, non-negative integer vector, length 6: `N_X`, `I_X`, `C_X`
#' for `X` in `h` (humans) and `d` (dogs), corresponding to the number of
#' susceptible (`N`), currently infected (`I`), and ever infected (`C`)
#' individuals in each population.
#'
#' There is always a single infected dog and no infected humans.
#'
#' @examples
#' dynfever_pop # to see function contents
#' dfpop <- dynfever_pop()
#' dfpop
#' dfpop <- dynfever_pop(dogs = 1e5)
#' dfpop
#' # note that non-sense arguments will produce errors
#' try(dynfever_pop(dogs = -1))
#' try(dynfever_pop(humans = 1.5))
#'
#' @family dynfever
#'
#' @export
dynfever_pop <- function(dogs = 1000L, humans = 1000L) {
  dogs |> check_integer() |> check_positive() |> check_scalar()
  humans |> check_integer() |> check_positive() |> check_scalar()
  c(N_d = dogs - 1L, I_d = 1L, C_d = 1L, N_h = humans, I_h = 0, C_h = 0)
}

#' @title Setup Dynamical Fever Parameters
#'
#' @param R0 a positive number: $R_0$, the reproductive number of Dynamical
#' Fever (DF) in a naive dog population
#'
#' @param pp a probability: the pairwise probability of dog-human effective
#' contact
#'
#' @param vax_humans a probability: the vaccinated fraction of the human
#' population
#'
#' @param vax_dogs a probability: the vaccinated fraction of the dog population
#'
#' @seealso [dynfever_solve()]
#'
#' @return a list of `numeric`s:
#'  - `R0`, the supplied $R_0$
#'  - `pp`, the supplied effective contact probability
#'  - `vax_dogs` and `vax_humans`, the vaccinated proportions
#'  - `p_avoid`, the probabilities for humans and dogs of avoiding DF
#'
#' @examples
#' dynfever_params # to see function contents
#' dfpar <- dynfever_params()
#' dfpar
#' dfpar <- dynfever_params(vax_dogs = 0.5, R0 = 3)
#' dfpar
#' # note, will complain about nonsensical parameters:
#' try(dynfever_params(R0 = -1))
#' try(dynfever_params(pp = 1.5))
#' try(dynfever_params(vax_humans = -1.5))
#' try(dynfever_params(vax_dogs = 1.5))
#'
#' @export
#' @family dynfever
dynfever_params <- function(
  R0 = 2, pp = 0.002, vax_humans = 0, vax_dogs = 0, N_d = dynfever_pop()["N_d"]
) {
  N_d |> check_integer() |> check_positive() |> check_scalar()
  R0 |> check_numeric() |> check_positive() |> check_scalar() |> check_lt(N_d)
  pp |> check_numeric() |> check_probability() |> check_scalar()
  vax_humans |> check_numeric() |> check_probability() |> check_scalar()
  vax_dogs |> check_numeric() |> check_probability() |> check_scalar()

  # transform transmission pars (R0, pp) => probabilities of avoiding infection
  p_avoid <- unname(c(1 - pp, 1 - R0 / N_d))

  return(as.list(environment()))
}

#' @title Single Dynamical Fever Step
#'
#' @description A function to compute the next state for the Dynamical Fever
#' Model; function signature to match the style of functions required by
#' [deSolve::ode()] family.
#'
#' @param t the time step
#'
#' @param y a named integer vector:
#' `c(N_human, N_dog, I_human, I_dog, C_human, C_dog)`
#'
#' @param parms a list: the model parameters; transformed from
#' [dynfever_parameters()] to `p_avoid = c(human = ..., dog = ...)` in
#' [dynfever_simulate()]
#'
#' @param ... ignored
#'
#' @details n.b. arguments are unchecked, as this is intended for internal use
#' in [dynfever_simulate()] repeatedly. The arguments are checked in that scope.
#'
#' @return a list: conforming to [deSolve::ode()] return values, a
#' `list(c(dN, dI, dC), FOI)`
#'
#' @examples
#' y <- dynfever_pop()
#' ps <- dynfever_params()
#' dy <- dynfever_dStep(y = y, parms = ps)
#' dy
#' y <- y + dy$states
#' dynfever_dStep(y = y, parms = ps)
#' # note, uses random number generator, so outcomes may differ:
#' dynfever_dStep(y = y, parms = ps)
#'
#' @export
#' @family dynfever
dynfever_dStep <- function(
  t, y, parms, ...
) with(parms, {
  dy <- y # copy y to carry over names / order; will replace all values
  # each time step, there is a binomial draw against susceptibles for exposure
  FOI <- setNames(1 - p_avoid^y["I_d"], c("FOI_h", "FOI_d"))
  # compute infection; note: this draws on the random number generator
  # (via runif) even e.g. if no susceptibles available
  infections <- qbinom(runif(c(1, 1)), size = y[c("N_h", "N_d")], prob = FOI)
  # infection removes from susceptibles
  dy[c("N_h", "N_d")] <- -infections
  # adds to cases (for cumulative counts)
  dy[c("C_h", "C_d")] <- infections
  # adds to infectious; but all currently infectious recover
  dy[c("I_h", "I_d")] <- infections - y[c("I_h", "I_d")]

  return(list(states = dy, params = FOI))
})

#' @title Simulate Dynamical Fever
#'
#' @description A function to simulate a Dynamical Fever (DF) outbreak; applies
#' [dynfever_dStep()] to iteratively generate the outbreak
#'
#' @param y a named vector of integers: the initial population of humans & dogs,
#' both infectious and susceptible
#'
#' @param times a series of integers, starting with 1L: the times to report;
#' result with also include time 0L.
#'
#' @param parms a list: the model parameters; see [dynfever_params()] return
#' value; the require entries are `vax_humans`, `vax_dogs`, `p_avoid`
#'
#' @param seed optionally, a random seed: if present, used with [set.seed()]
#'
#' @details This method is not exported, but can still be viewed by students via
#' `ICI3D:::dynfever_simulate` (note the *triple* `:::`)
#'
#' @return a `data.table` with columns `time`, ... TODO
#'
#' @importFrom data.table rbindlist
#'
#' @examples
#' library(data.table) # for pretty-printing output
#' set.seed(1); ps <- dynfever_params(); pop <- dynfever_pop()
#' sim1 <- ICI3D:::dynfever_simulate(pop, 1:10, ps)
#' sim1$states
#' sim1$params
#' # note that dynfever_simulate() is stochastic:
#' sim2 <- ICI3D:::dynfever_simulate(pop, 1:10, ps)
#' !all(sim1$states == sim2$states)
#' # however, it can be provided a seed to ensure comparable results:
#' all(sim1$states == ICI3D:::dynfever_simulate(pop, 1:10, ps, 1L)$states)
#'
#' @export
#' @family dynfever
dynfever_simulate <- function(y, ts, ps, seed) {
  if (!missing(seed)) set.seed(seed)

  # apply vaccination; note: this draws on the random number generator
  # (via runif) even e.g. if vaccinated proportions are 0.
  vs <- with(ps, { # draw vaccinated count
    qbinom(runif(c(1, 1)), size = y[c("N_h", "N_d")], prob = c(vax_humans, vax_dogs))
  }) |> setNames(c("V_h", "V_d"))
  # remove vaccinees from susceptibles
  y[c("N_h", "N_d")] <- y[c("N_h", "N_d")] - vs

  # setup initial population / parameters
  yt <- list(
    states = y,
    params = c(NA_real_, NA_real_) |> setNames(c("FOI_h", "FOI_d"))
  )

  # setup the storage `data.table`s for "states" and "params"
  states <- within(data.table(time = c(0L, ts), key = "time"), {
    N_h <- V_h <- I_h <- C_h <- NA_integer_
    N_d <- V_d <- I_d <- C_d <- NA_integer_
  })
  # row 1 => time == 0
  states[1, names(yt$states) := as.list(yt$states) ]
  states[, names(vs) := as.list(vs) ]
  params <- within(data.table(time = c(0L, ts), key = "time"), {
    FOI_h <- FOI_d <- NA_real_
  })

  ti <- 1L; has_FOI <- TRUE
  # while there is FOI AND still time remaining ...
  while (has_FOI && (ti <= length(ts))) {
    # get the dy / parameters
    nexty <- dynfever_dStep(ti, yt$states, ps)
    # update the loop variables & record results
    yt$states <- yt$states + nexty$states; yt$params <- nexty$params
    ti <- ti + 1L; has_FOI <- sum(yt$params) != 0
    states[ti, names(yt$states) := as.list(yt$states) ]
    params[ti, names(yt$params) := as.list(yt$params) ]
  }

  # fill in any leftover outcomes by carrying forward end of epidemic
  setnafill(states, type = "locf")
  setnafill(params, type = "locf")

  return(list(states = states, params = params))
}

#' @title Dynamical Fever Sampler
#'
#' @param n a positive integer: the number of samples
#'
#' @param ref.seed see [set.seed()]; the initial seed. each sample will be a
#' specific shift from this seed
#'
#' @param parms a list: either same as return of [dynfever_params()] or a list
#' of elements that are each the same as return of [dynfever_params()]
#'
#' @inheritParams dynfever_simulate
#'
#' @seealso [dynfever_simulate()]
#'
#' @return a  pair of `data.table`s, same as [dynfever_simulate()] with added
#' a `sample` columns.
#'
#' @examples
#' one_run <- dynfever_simulate(seed = 1L)
#' ten_run <- dynfever_sample(n = 10, ref.seed = 1L)
#' # note:
#' all(one_run$states == ten_run$states[sample == 1][, -1])
#' # dynfever_sample() can also be used with a series of paramters
#' pars <- seq(.1, .5, by = .1) |> lapply(
#'  \(vax_d) dynfever_params(vax_dogs = vax_d)
#' )
#'
#' five_run <- dynfever_sample(n=5, parms = pars)
#'
#' @export
#' @importFrom parallel mcmapply
#' @family dynfever
dynfever_sample <- function(
  n = 1L, y = dynfever_pop(),
  times = 1L:MAXTIME,
  parms = dynfever_params(N_d = y["N_d"]),
  MAXTIME = 30, ref.seed,
  ...
) {
  times |> check_integer() |> check_series() |> check_positive()
  y |> check_integer() |> check_nonnegative()
  names(y) |> check_among(c("N_h", "I_h", "C_h", "N_d", "I_d", "C_d"))
  # note, no check of `n` required, as seq_seeds / seq_parms does that check

  # size seed and params as requested
  seeds <- seq_seeds(ref.seed, n)
  parms <- seq_parms(parms, n)

  smpls <- parallel::mcmapply(
    dynfever_simulate,
    ps = parms, seed = seeds, MoreArgs = list(y = y, ts = times),
    SIMPLIFY = FALSE
  )

  refkeys <- smpls[[1]] |> lapply(\(dt) c("sample", key(dt)))

  res <- names(refkeys) |> lapply(
    \(part) smpls |> lapply(\(r) r[[part]]) |> rbindlist(idcol = "sample") |> setkeyv(refkeys[[part]]) |> setcolorder()
  ) |> setNames(names(refkeys))

  return(res)

}

#' @title Long-Format Dynamical Fever Output
#'
#' @description Reformat [dynfever_simulate()] and [dynfever_sample()] output
#' to "long" format.
#'
#' @param dynfever_output the output of [dynfever_simulate()] and [dynfever_sample()];
#' either the full list, with elements `states` and `params`, both
#' `data.table`s; or directly one of those two elements. If the input is already
#' in long format, it will be returned unchanged. See [dynfever_simulate()] and
#' [dynfever_sample()] for the base output formats
#'
#' @return depending on input:
#'  - if given a list, return a list with elements `states` and `params`, both
#' `data.table`s,
#'  - if given a `data.table`, return a `data.table`
#'
#' In both cases, the `data.table`(s) have the keys of their input: optionally
#' `sample` (when given output from [dynfever_sample()], always `time` and
#' `pop`, then `state` (for `states` output) or `param` (for `params` output)
#'
#' @examples
#' require(data.table)
#' dyn_one <- dynfever_simulate()
#' dyn_run <- dynfever_sample(n=10)
#'
#' one_long <- dynfever_toLong(dyn_one)
#' one_long
#' # also works for just e.g. states:
#' one_long_states <- dynfever_toLong(dyn_one$states)
#' all(one_long$states == one_long_states)
#' # and for multiple samples:
#' run_long <- dynfever_toLong(dyn_one)
#' run_long
#' # note preservation of keys:
#' key(one_long$states)
#' key(run_long$states)
#'
#' @export
#' @family dynfever
dynfever_toLong <- function(dynfever_output) {

  # if we have the pair of outputs ...
  if (!is.data.table(dynfever_output)) { # apply to both
    return(dynfever_output |> lapply(dynfever_toLong))
  } else { # otherwise ...
    skey <- key(dynfever_output)
    # if already wide, do nothing ...
    if (length(setdiff(skey, c("sample", "time"))) != 0) {
      return(dynfever_output)
    } else { # otherwise, determine if states vs params
      target <- if (any(names(dynfever_output) %like% "N_d")) "state" else "param"
      # convert to long format
      long <- melt.data.table(
        dynfever_output, id.vars = skey,
        value.name = if (target == "state") "count" else "value"
      )
      # split variable into target and pop, then remove it
      long[, c(target, "pop") := tstrsplit(variable, "_")]
      long$variable <- NULL
      # reformat the pop variable
      long[, pop := factor(pop, levels = c("d", "h"), ordered = TRUE) ]
      # reformat the target variable
      if (target == "state") {
        long[, state := factor(state, levels = c("N", "V", "I", "C"), ordered = TRUE) ]
      } else {
        long[, param := factor(param, ordered = TRUE) ]
      }
      return(long |> setkeyv(c(skey, "pop", target)) |> setcolorder())
    }
  }

}

#' @title Summarize Dynamical Fever Results
#'
#' @description Reports the peak value and timing for infections, as well as the
#' final size and outbreak duration, stratified by population.
#'
#' @inheritParams dynfever_toLong
#'
#' @export
#' @family dynfever
dynfever_summarize <- function(dynfever_output) {
  long.dts <- dynfever_toLong(dynfever_output)
  long.dts$summary <- long.dts$states[state %in% c("I", "C"), .SD[which.max(count)], by=setdiff(key(long.dts$states), "time")]
  return(long.dts)
}



