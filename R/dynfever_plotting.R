
#' @title Convert Dynamical Fever Output to Incidence
#'
#' @description Takes raw [dynfever_simulate()] and [dynfever_sample()] output
#' and extracts incidence
#'
#' @inheritParams dynfever_toLong
#'
#' @export
#' @family dynfever
dynfever_incidence <- function(dynfever_output) {
  # if necessary, get the relevant part of the outputS
  if (!is.data.table(dynfever_output)) {
    states.dt <- dynfever_output$states
  } else {
    states.dt <- dynfever_output
  }
  names(states.dt) |> check_contains(c("time", "N_d"))

  # convert to long format
  states.dt <- states.dt |> dynfever_toLong()

  # only interested in the cumulative infections => their `diff`s
  res <- states.dt[state == "C"][, .(
      time, incidence = c(count[1], diff(count))
    ),
    by = setdiff(key(states.dt), "time") # group by the normal keys, except time
  ] |> setkeyv(key(states.dt)) # re-establish keys

  return(res)
}

#' @title Dynamical Fever Dog & Human Incidence Plot
#'
#' @param sample_id an integer: a particular sample from outputs to consider
#'
#' @inheritParams dynfever_toLong
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @export
#' @family dynfever
dynfever_plot_simple <- function(
  dynfever_output = dynfever_sample(), sample_id
) {
  # extract plotting data & prepare filter
  summary.dt <- dynfever_summarize(dynfever_output)$summary
  incidence.dt <- dynfever_incidence(dynfever_output)
  if (!missing(sample_id)) {
    filterfun = \(dt) subset(dt, sample == sample_id)
  } else {
    filterfun = \(dt) dt
  }
  # create plot
  ggplot(incidence.dt) + aes(x = time + 1 - 0.5, y = incidence) +
    facet_grid( # facet by populations ...
      cols = vars(pop), labeller = labeller(pop = c(d="Dogs", h="Humans"))
    ) +
    geom_label(
      aes(label = label),
      data = \(dt) data_grid_labels(dt, incidence, time = 2),
      color = "grey20", fill = alpha("white", .5), label.size = 0
    ) +
    geom_col(width = 1, fill = alpha("darkgrey", .5), color = "black", data = filterfun) +
    geom_label_repel(aes(
      label = sprintf(
        "%s: %i%s week %i",
        fifelse(state == "I", "Peak", "Total"), count, fifelse(state == "I", " in", "\nby"), time
      )
    ), data = filterfun(summary.dt)[,.(
      pop, incidence = fifelse(state == "I", count, 0)+3,
      state, count, time
    )]) +
    scale_y_continuous("Cases") +
    scale_x_continuous(NULL) +
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_minimal() + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text = element_text(size = rel(1))
    )
}

#' @rdname dynfever_plot_simple
#' @export
#' @family dynfever
dynfever_plot_empty <- function() {
  # create plot
  ggplot(dynfever_incidence(dynfever_DAIDD_outbreak)) +
    aes(x = time + 1 - 0.5, y = incidence) +
    facet_grid( # facet by populations ...
      cols = vars(pop), labeller = labeller(pop = c(d="Dogs", h="Humans"))
    ) +
    geom_label(
      aes(label = label),
      data = \(dt) data_grid_labels(dt, incidence, time = 2),
      color = "grey20", fill = alpha("white", .5), label.size = 0
    ) +
    geom_blank() +
    scale_y_continuous("Cases") +
    scale_x_continuous(NULL) +
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_minimal() + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      strip.text = element_text(size = rel(1))
    )
}

#' @title Plot DynFever Outputs as an Annual Series
#'
#' @description Plots a multi-year Dynamical Fever series
#'
#' @param dynfever_output see return value of [dynfever_sample()]
#'
#' @param startyear the "year 0" corresponding to `sample == 1` outcomes
#'
#' @export
#' @family dynfever
dynfever_plot_series <- function(
  dynfever_output = dynfever_sample(n = 5),
  startyear = {
    res <- Sys.Date() |> as.POSIXlt()
    res$year <- res$year - 5
    res$mon <- 0
    res$mday <- 1
    res
  }
) {
  if (!is(startyear, "POSIXlt")) {
    startyear |> check_integer()
    startyear <- startyear |> ISOdate(year = _, month = 1, day = 1) |> as.POSIXlt()
  }
  # get the incidence from the outputs
  incidence.dt <- dynfever_incidence(dynfever_output)
  incidence.dt[, date := {
    bs <- startyear
    # sample -> relative year to start date (though indexed from 1)
    bs$year <- bs$year + sample - 1
    # time in weeks, 0 indexed
    (bs |> as.Date()) + time * 7
  } ]

  ggplot(incidence.dt) + aes(x = date, y = incidence) +
    facet_grid( # facet by populations ...
      rows = vars(pop), labeller = labeller(pop = c(d="Dogs", h="Humans")),
      switch = "y"
    ) +
    geom_label(
      aes(label = label),
      data = \(dt) data_grid_labels(dt, incidence, date = startyear |> as.Date() + 26*7),
      color = "grey40", fill = alpha("white", .5), label.size = 0
    ) +
    geom_col(width = 7, fill = "darkgrey") +
    scale_y_continuous("Cases in ...") +
    scale_x_date(NULL, date_breaks = "years", date_label = "%Y") +
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_minimal() + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text = element_text(size = rel(1)),
      strip.placement = "inside", panel.spacing.y = unit(1, "line"),
      axis.text.x = element_text(hjust = 0, size = rel(1.25), color = "black")
    )

}

#' @title Dynamical Fever: Compare Vaccine Pops
#'
#' @param dt a [data.table::data.table()]; with columns:
#'  - `sample`, the sample id
#'  - `vax_pop`, the vaccinated population, either "h" or "d"
#'  - `pop`, the measured population, either "h" or "d"
#'  - `vax_cov`, the coverage in the vaccinated population
#'  - `count`, the final size of the outbreak in `pop`
#'
#' @param countbinwidth positive integer, how wide to make the bins for
#' counts?
#'
#' @return a [ggplot2::ggplot()] object
#' @export
#' @family dynfever
dynfever_plot_heatmap <- function(
  dt, countbinwidth = 10
) {

  # bin outcomes, then get their relative frequencies
  heat_dt <- dt[,
    .N, keyby = .(
      vax_cov, vax_pop, pop,
      count = floor(count/countbinwidth)*countbinwidth + countbinwidth/2
    )
  ][, .(
      count, frac = N/sum(N)
    ), keyby = .(vax_cov, vax_pop, pop)
  ]

  # get the median outcomes
  med_dt <- dt[, .(
    count = median(count)
  ), by = .(vax_cov, vax_pop, pop)]

  # get center point references
  center_dt <- dt[,
    .(count = mean(range(count))), by = pop
  ][, dt[,
      .(vax_cov = mean(range(vax_cov))), by = vax_pop
  ], by = .(pop, count)][, needed := TRUE ]

  # prepare to display "???" for missing ones
  center_dt[
    heat_dt[, .N, by=.(pop, vax_pop)], on = .(pop, vax_pop),
    needed := !(N > 0)
  ]

  return(ggplot(heat_dt) +
    aes(vax_cov, count, fill = log(frac)) + facet_grid(
      pop ~ vax_pop, labeller = labeller(
        pop = c(h = "People", d = "Dogs"),
        vax_pop = c(h = "People", d = "Dogs")
      ),
      switch = "y"
    ) +
    geom_tile(alpha = 0.7) +
    geom_point(
      aes(fill = NULL, color = "median"),
      data = med_dt
    ) + coord_cartesian(clip = "off") +
    geom_text(
      aes(fill = NULL),
      data = center_dt[needed == TRUE],
      label = "???", size = 10
    ) +
    scale_x_continuous("Coverage in ...", expand = c(0, 0), position = "top") +
    scale_y_continuous("Incidence in ...", expand = c(0, 0)) +
    scale_fill_continuous(
      "Fraction\nof Simulations",
      labels = scales::math_format(10^.x)
    ) +
    scale_color_manual(
      NULL, labels = c(median = "Median"), values = c(median = "red")
    ) +
    theme_minimal() + theme(
      panel.spacing = unit(1.5, "line"), strip.placement = "outside",
      legend.position = "bottom"
    ))
}
