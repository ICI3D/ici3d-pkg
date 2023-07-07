#' @import data.table ggplot2 shiny
NULL

#' @title Heterogeneous Population Generation
#'
#' @description Generates a heterogeneous population of contact rates
#'
#' @param n an integer: the number of individuals to generate
#'
#' @param beta_mean a numeric, greater than 0: mean contact rate
#'
#' @param beta_var a numeric, greater than 0: variance of the contact rate
#'
#' @param seed an integer, optional: the random seed to use (if supplied); see
#' [set.seed()]
#'
#' @return a numeric vector, length `n`: contact rates, sorted in descending
#' order
#'
#' @details This function generates a heterogeneous population of contact rates
#' using [rgamma()].  The mean and variance of the distribution are
#' transformed into the `shape` and `scale` parameters for [rgamma()].
#'
#' @family heterogeneity
#'
#' @examples
#' mean_contacts <- 2
#' hp <- het_population(n = 1000, beta_mean = mean_contacts, beta_var = .1)
#' het_contact_hist(hp, beta_mean = mean_contacts)
#'
#' @export
het_population <- function(
  n, beta_mean, beta_var, seed
) {
  # ensure arguments are appropriate
  n |> check_integer() |> check_positive()
  beta_mean |> check_numeric() |> check_positive()
  beta_var |> check_numeric() |> check_positive()

  # set the seed if desired
  if (!missing(seed)) set.seed(seed)

  # transform the arguments, draw gamma deviates, sort from high => low
  return(
    rgamma(n, shape = beta_mean^2 / beta_var, scale = beta_var / beta_mean) |>
    sort(decreasing = TRUE)
  )
}

#' @title Heterogeneous Contact Rate Histogram
#'
#' @description Plots a histogram of contact rates
#'
#' @param contact_rates numeric vector; contact rates to plot, e.g. from
#' [het_population()]
#'
#' @inheritParams het_population
#' @inheritParams ggplot2::geom_histogram
#' @inheritDotParams ggplot2::geom_histogram
#'
#' @return a ggplot object
#' @family heterogeneity
#' @examples
#' mean_contacts <- 2
#' hp <- het_population(n = 1000, beta_mean = mean_contacts, beta_var = .1)
#' het_contact_hist(hp, beta_mean = mean_contacts)
#'
#' @export
het_contact_hist <- function(
  contact_rates = numeric(), beta_mean, binwidth = .1, ...
) {
  # check arguments
  contact_rates |> check_numeric() |> check_positive()
  beta_mean |> check_numeric() |> check_positive() |> check_scalar()

  # create a minimum scale limit
  pseudomax <- 2*beta_mean

  # create a plot ...
  return(
    ggplot(data.table(contact_rates = contact_rates)) + # from contact rates...
    geom_blank( # at least shown on c(0, pseudomax) ...
      aes(lims),
      data.table(lims = c(0, pseudomax))
    ) +
    geom_histogram( # showing a histogram of the contact rates ...
      aes(x = contact_rates, y = after_stat(density)),
      binwidth = binwidth, ...
    ) +
    geom_vline(xintercept = beta_mean, color = "red") + # indicating the mean...
    scale_y_continuous(name = NULL, guide = "none") + # don't show y axis
    scale_x_continuous(
      name = "contact rate (1/day)",
      expand = expansion()
    ) +
    theme_bw() + theme( # minimal theme + bigger text
      panel.border = element_blank(),
      axis.title = element_text(size = rel(2)),
      axis.text = element_text(size = rel(2))
    )
  )
}

#' @title Heterogeneous Final Size Histogram
#'
#' @description Plots a histogram of final sizes from many [het_sample()]s
#'
#' @param ts a [data.table::data.table()]; the return from [het_sample()]
#'
#' @return a [ggplot2::ggplot()] object
#' @family heterogeneity
#' @examples
#' hp <- het_population(n = 500, beta_mean = 2, beta_var = 1)
#' samples <- het_sample(100, hp)
#' het_finalsize_hist(samples, 10)
#'
#' @export
het_finalsize_hist <- function(ts, binwidth = 5, ...) {
  maxpop <- ts[1, S+I+R] # get the max total population
  # extract the final sizes; note, since het_simulate allows immune waning, have
  # to use accumulator state rather than the final `R` value
  plot_ts <- ts[, .(value = cumulativeI[.N]), by = sample_id ]

  # plot a histogram (y-axis categories, x-axis counts)
  return(
    ggplot(plot_ts) +
    aes(y = value) + geom_histogram(binwidth = binwidth, ...) +
    theme_minimal() + theme(
      axis.title = element_text(size=rel(2)),
      axis.text = element_text(size=rel(2))
    ) + coord_cartesian(ylim = c(0, maxpop), expand = FALSE, clip = "off") +
    labs(y = "... of size", x = "# of outbreaks")
  )
}

#' @title Simulate Heterogeneous Transmission Model
#'
#' @description
#' Run the heterogeneous transmission model for a particular parameter
#' combination.
#'
#' @param mxdst individual contact rates, as returned by [het_population()]
#'
#' @param params a list, with `tmax`, `recovery_rate`, and `waning_rate` keys
#'
#' @importFrom data.table between
#' @family heterogeneity
#' @export
het_simulate <- function(
  mxdst, params, seed
) {
  if (!missing(seed)) { set.seed(seed) }
  with(params, {

  pop.size <- length(mxdst)

  # create individual state trackers: clocks for recovery & waning
  trec <- rep(Inf, pop.size)
  twan <- -trec
  # ... and who's infectious
  whichinf <- rep(FALSE, pop.size)

  # an individual is susceptible if:
  # !whichinf & (twan < 0)

  # sample initial infectee, set recovery + waning times
  tee <- sample.int(pop.size, 1, prob = mxdst)
  whichinf[tee] <- TRUE
  trec[tee] <- rexp(1, recovery_rate)

  if (exists("waning_rate") && (waning_rate > 0)) {
    wan_sample <- function() rexp(1, waning_rate)
  } else {
    wan_sample <- function() Inf
  }

  twan[tee] <- trec[tee] + wan_sample()

  # create loop variables
  con.rate <- sum(mxdst[whichinf])
  icur <- 1
  tcur <- 0
  len <- 1

  # create containers to update with changes when events occur; these are Deltas
  event.times <- list(0)
  S <- list(pop.size-1)
  I <- list(1)
  R <- list(0)
  cumulativeI <- list(1)

  while ((icur > 0) & (tcur < tmax)) {
    nextt <- rexp(1, con.rate)
    isrecovery <- any(nextt > trec)

    # if a recovery would happen first, amend event time to recovery ...
    if (isrecovery) {

      # find the individual + associated recovery time ...
      irec <- which.min(trec)
      nextt <- trec[irec]

      # update the infectious individuals + recompute contacting rate
      whichinf[irec] <- FALSE
      trec[irec] <- Inf
      con.rate <- sum(mxdst[whichinf])

      # update loop control information
      icur <- icur - 1

    }

    # update all the timers accordingly, and identify any waned protection
    tcur <- tcur + nextt
    trec <- trec - nextt
    wanes <- between(twan, 0, nextt)
    want <- sort(twan[wanes]) # for aggregation, we care about timing, not id
    twan <- twan - nextt


    # wanes => R -> S; if any waning:
    if (length(want)) {

      # change to relative delays
      nextt <- nextt - want[length(want)] # happens some Delta after last want
      want <- c(want[1], diff(want))

      # make the new indices + update len
      app <- len + seq_along(want)
      len <- len + length(want)

      # update Delta containers; n.b. [] vs [[]] access for slicing
      R[app] <- -1 # move from R
      S[app] <- 1  # to S
      event.times[app] <- want # at relative times
      # no changes (yet) to I
      I[app] <- 0
      cumulativeI[app] <- 0

    }

    if (isrecovery) {
      # recovery an I to R at the relative time; no change to cumulative inc, S
      I[[len + 1]] <- -1
      R[[len + 1]] <- 1
      event.times[[len + 1]] <- nextt

      S[[len + 1]] <- 0
      cumulativeI[[len + 1]] <- 0
      len <- len + 1
    } else {
      # choose infector
      tor <- which(whichinf)[sample(icur, size = 1, prob = mxdst[whichinf])]
      # choose infectee - draw from n-1, then shift if necessary
      tee <- sample(pop.size - 1, 1, prob = mxdst[-tor])
      if (tee >= tor) tee <- tee + 1

      # if tee infectable => infect
      if (!whichinf[tee] & (twan[tee] < 0)) {
        # infect tee, update the infectious contact rate
        whichinf[tee] <- TRUE
        trec[tee] <- rexp(1, rate = recovery_rate)
        twan[tee] <- trec[tee] + wan_sample()
        con.rate <- sum(mxdst[whichinf])

        # update state records
        S[[len + 1]] <- -1
        I[[len + 1]] <- 1
        R[[len + 1]] <- 0
        event.times[[len + 1]] <- nextt
        cumulativeI[[len + 1]] <- 1
        len <- len + 1

        # update loop variable
        icur <- icur + 1
      } # otherwise, no change
    }
    # done with one time draw ...
  } # end simulation loop

  # potentially extend to tmax ...
  day <- cumsum(event.times)
  if (day[length(day)] < tmax) {
    day <- c(day, tmax)
    S[[len+1]] <- 0
    I[[len+1]] <- 0
    R[[len+1]] <- 0
    cumulativeI[[len+1]] <- 0
  }

  # convert deltas to time series, filter to tmax
  return(data.table(
    day = day,
    S = cumsum(S), I = cumsum(I), R = cumsum(R),
    cumulativeI = cumsum(cumulativeI)
  )[ day <= tmax ])
})
}

#' @title Time Series Plot for Heterogeneity Tutorial
#'
#' @description
#' Produces a standard compartment time series plot for outputs from
#' [het_sample()].
#'
#' @param ts a [data.table()], as yielded by [het_sample()]
#'
#' @return a [ggplot2::ggplot()] object
#'
#' @export
#' @family heterogeneity
#' @examples
#' het_plot(het_sample(30))
#'
het_plot <- function(ts) {
  tmax <- ts[, max(day)]
  maxpop <- ts[1, S + I + R]
  n <- ts[, max(sample_id)]
  long_dt <- ts[, .SD, .SDcols = -patterns("^cumulative")] |>
    data.table::melt.data.table(
      id.vars = c("sample_id", "day"), variable.name = "state"
    )
  rug_dt <- ts[,
    if (any(I == 0)) .SD[I==0][1] else .SD[.N], by = sample_id
  ][, .SD, .SDcols = -patterns("^cumulative")] |>
    data.table::melt.data.table(
      id.vars = c("sample_id", "day"), variable.name = "state"
    )
  return(ggplot(long_dt) +
    aes(
      x = day, y = value,
      color = state, group = interaction(sample_id, state)
    ) +
    geom_step(alpha = sqrt(1/n)) +
    geom_rug(
      data = rug_dt,
      alpha = (1/n)^(1/3), outside = TRUE,
      sides = "r"
    ) +
    geom_rug(
      data = subset(rug_dt, state == "I"),
      alpha = (1/n)^(1/3), outside = TRUE,
      sides = "t"
    ) +
    coord_cartesian(
      xlim = c(0, tmax),
      ylim = c(0, maxpop),
      clip = "off",
      expand = FALSE
    ) +
    scale_x_continuous("Day") +
    scale_y_continuous("Individuals") +
    scale_color_discrete(
      name = NULL, guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    scale_alpha_continuous(guide = "none", range = c(0.05, 1)) +
    theme_minimal() + theme(
      axis.title = element_text(size = rel(2)),
      axis.text = element_text(size = rel(2)),
      plot.margin = margin(t = 1, r = 1, unit = "line")
    ) +
    theme(legend.position = "left", legend.direction = "vertical")
  )
}

#' @title Sample [het_simulate()]
#'
#' @description
#' Run [het_simulate()] many times with the same parameter combinations.
#'
#' @param n a positive integer, the number of samples to simulate
#'
#' @inheritParams het_simulate
#'
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist melt.data.table
#' @family heterogeneity
#' @examples
#' hetruns <- het_sample(10)
#' hetruns[day == 5.0]
#'
#' @export
het_sample <- function(
  n,
  mxdst = het_population(n = 300, beta_mean = 1, beta_var = .5),
  params = list(tmax = 5, recovery_rate = 1, waning_rate = 0)
) {

  # check input parameters
  n |> check_integer() |> check_positive() |> check_scalar()
  mxdst |> check_numeric() |> check_positive()

  # create a sample frame => for each sample, run het_simulate
  ts <- seq_len(n) |> parallel::mclapply(
    FUN = \(sid, mxdst, params) het_simulate(mxdst, params, sid),
    mxdst = mxdst, params = params,
    mc.preschedule = FALSE
  ) |> rbindlist(idcol = "sample_id")

  return(ts |> setkey(sample_id, day))

}
