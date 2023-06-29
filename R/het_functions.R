#' @import data.table ggplot2 shiny
NULL

#' @title Heterogeneous Population Generation
#'
#' @description Generates a heterogeneous population of contact rates
#'
#' @param n an integer: the number of individuals to generate
#'
#' @param beta_mean a numeric: the mean contact rate
#'
#' @param beta_var a numeric: the variance of the contact rate
#'
#' @param seed an integer, optional: the random seed to use (if supplied)
#'
#' @return a numeric vector of length `n` of contact rates, sorted in descending
#' order
#'
#' @details This function generates a heterogeneous population of contact rates
#' using a gamma distribution.  The mean and variance of the distribution are
#' transformed into the shape and scale parameters for [rgamma()].
#'
#' @family heterogeneity
#'
#' @examples
#' hp <- het_population(n = 1000, beta_mean = 1, beta_var = 1)
#' hist(hp)
#'
#' @export
het_population <- function(
  n, beta_mean, beta_var, seed
) {
  if (!missing(seed)) set.seed(seed)
  scale <- beta_var / beta_mean
  shape <- beta_mean / scale
  rgamma(n, shape = shape, scale = scale) |> sort(decreasing = TRUE)
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
#' hp <- het_population(n = 1000, beta_mean = 1, beta_var = 1)
#' het_plot_hist(hp, beta_mean = 1)
#'
#' @export
het_contact_hist <- function(
  contact_rates, beta_mean, binwidth = .1, ...
) {
  return(ggplot(data.table(contact_rates = contact_rates)) +
    geom_histogram(
      aes(x = contact_rates, y = after_stat(density)),
      binwidth = binwidth, ...
    ) + xlab("contact rate (1/day)") +
    geom_vline(xintercept = beta_mean, color = "red") +
    theme_bw() + theme(
      panel.border = element_blank(),
      axis.title = element_text(size = rel(2)),
      axis.text = element_text(size = rel(2))
  ))
}

#' @title Heterogeneous Final Size Histogram
#'
#' @description Plots a histogram of final sizes from many [het_sample()]s
#'
#' @param dt a data.table; the return from [het_sample()]
#'
#' @return a ggplot object
#' @family heterogeneity
#' @examples
#' hp <- het_population(n = 1000, beta_mean = 1, beta_var = 1)
#' samples <- het_sample(hp)
#' het_finalsize_hist(samples)
#'
#' @export
het_finalsize_hist <- function(ts, binwidth = 5, ...) {
  maxpop <- ts[1, S+I+R]
  return(
    ggplot(
      ts[, .(value = cumulativeI[.N]), by = sample_id ]
    ) + aes(x = value) + geom_histogram(binwidth = binwidth) +
    theme_minimal() + theme(
      axis.title = element_text(size=rel(2)),
      axis.text = element_text(size=rel(2))
    ) + coord_cartesian(xlim = c(0, maxpop)) +
    labs(x= "... of size", y = "# of outbreaks")
  )
}

het.runs.hist <- function(ts) ggplot(
  ts[,list(`final size` = max(cumulativeI)), by=runid]
) + aes(x=`final size`) + geom_histogram(binwidth = 5) +
  theme_minimal() + theme(
    axis.title = element_text(size=rel(2)),
    axis.text = element_text(size=rel(2))
  ) + ylab("# of outbreaks")

#' @title Sample Pairs of Distinct Indices
#'
#' @param n positive integer: number of pairs to sample
#'
#' @param prob numeric vector: the population weights for sampling each index;
#' if length(prob) == 1, then treated as the population size and all indices
#' equally weighted
#'
#' @param verbose Print progress messages?
#'
#' @return A list with two elements, `source` and `sink`, each of length `n`,
#' corresponding to the sampled pairs indices
#'
#' @family heterogeneity
#' @examples
#' draws <- sample_pair(n = 100, prob = 10)
#'
#' @export
sample_pair <- function(n, prob, verbose = FALSE) {
  if (length(prob) == 1) {
    sampler <- function(sn) {
      sample(prob, sn, replace = TRUE)
    }
  } else {
    sampler <- function(sn) {
      sample(length(prob), sn, replace = TRUE, prob = prob)
    }
  }
  source <- sampler(n)
  sink <-   sampler(n)
  redraw <- which(source == sink)
  while (length(redraw)) {
    if (verbose) {
      message(sprintf("sample_pair: redrawing %i ...", length(redraw)))
    }
    sink[redraw] <- sampler(length(redraw))
    redraw <- redraw[which(source[redraw] == sink[redraw])]
  }
  return(list(source = source, sink = sink))
}

#' @title Sample Event Times Until a Maximum Time
#'
#' @param maxval a positive numeric: the maximum time to simulate until
#'
#' @param meantime a positive numeric: the mean time between events
#'
#' @param sampler a function: the sampling distribution for time between events
#'
#' @param ... any parameter arguments to `sampler`
#'
#' @param oversample a numeric: the amount to oversample by when generating
#'
#' @param verbose Print progress messages?
#'
#' @return A numeric vector of event times: the cumulative sum of draws from
#' `sampler` until the last value exceeds `maxval`. That sum is then trimmed to
#' values `<= maxval`. The first value in the vector is always 0.
#'
#' @family heterogeneity
#' @examples
#' draws <- sample_until(n = 100, meantime = 1, rexp, rate = 1)
#' cumsum(draws)
#'
#' @export
sample_until <- function(
  maxval, meantime, sampler, ...,
  oversample = 0.1, verbose = FALSE
) {
  # generate a vector of interaction times - oversample by 10%
  expectedevents <- maxval / meantime
  times <- c(
    0, cumsum(sampler(ceiling(expectedevents * (1 + oversample)), ...))
  )
  # make sure we have enough events - extend by 10% until we do
  while (times[length(times)] < maxval) {
    if (verbose) {
      message(sprintf(
        "sample_until: current end %d vs maxval %d -> add %i events ...",
        times[length(times)], maxval, ceiling(expectedevents * oversample)
      ))
    }
    times <- c(
      times,
      times[length(times)] +
        cumsum(sampler(ceiling(expectedevents * oversample), ...))
    )
  }
  # trim off any events that occur after maxval
  times <- times[times <= maxval]

  return(times)
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
#' @family heterogeneity
#' @export
het_simulate <- function(
  mxdst, params, seed
) {
  if (!missing(seed)) { set.seed(seed) }
  pop.size <- length(mxdst)
  end.time <- params$tmax
  gmma <- params$recovery_rate
  rho <- params$waning_rate

  tloss <- trec <- rep(NaN, pop.size)

  inf.rate <- sum(mxdst)

  event.times <- c(0, rexp(ceiling(end.time/inf.rate*1.1), inf.rate))
  while(sum(event.times) < end.time) event.times <- c(event.times, rexp(ceiling(end.time/inf.rate*0.1), inf.rate))

  event.times <- cumsum(event.times)
  event.times <- event.times[which(event.times <= end.time)]

  stepmax <- length(event.times)
  S <- c(pop.size-1, rep.int(0, stepmax-1))
  I <- c(1, rep.int(0, stepmax-1))
  R <- rep.int(0, stepmax)
  cumulativeI <- I

  new.ted <- sample.int(pop.size, 1, prob = mxdst)
  last.time <- 0 # first individual is infected at time 0

  ( last.time + rexp(1, gmma) -> trec[new.ted] ) +
    ifelse(rho, rexp(1, rho), end.time) -> tloss[new.ted]
  ## initialize time series data frame
  #ts <- data.table(day=new.time, S=pop.size-1, I=1, R = 0, cumulativeI=1)
  #next.time <- data.table::copy(ts)
  ## calculate current inf rate - should this be changing?  was inside while loop


  ## expected samples to get sum(samples) > end.time
  # times <- rexp(ceiling(end.time*1.1*inf.rate), inf.rate)
  # while (sum(times) < end.time) times <- c(times, rexp(ceiling(end.time*0.1*inf.rate), inf.rate))

  last.time <- 0
  for (i in 2:stepmax) {
    recoveries <- which(trec < last.time)
    trec[recoveries] <- NaN
    I[i] <- I[i-1] - length(recoveries)
    R[i] <- R[i-1] + length(recoveries)
    S[i] <- S[i-1]
    cumulativeI[i] <- cumulativeI[i-1]

    who.loss <- which(tloss < last.time)
    if (ilosses <- length(who.loss)) {
      tloss[who.loss] <- NaN
      R[i] <- R[i] - ilosses
      S[i] <- S[i] + ilosses
    }

    cur.tor <- sample(pop.size, 1, prob = mxdst)
    ## has this id been infected?
    if(!(is.nan(tloss[cur.tor]) | is.nan(trec[cur.tor]))) {
      # choose from everyone except infec*tor*
      new.ted <- sample((1:pop.size)[-cur.tor], 1, prob = mxdst[-cur.tor])
      if(is.nan(tloss[new.ted])) { ## is this ID currently susceptible?
        ## then infect them
        S[i] <- S[i] - 1
        I[i] <- I[i] + 1
        cumulativeI[i] <- cumulativeI[i] + 1
        # infected @ new.time, recover @ infected + recovery draw, loss @ recover + loss draw
        ( event.times[i] + rexp(1, gmma) -> trec[new.ted] ) +
          ifelse(rho, rexp(1, rho), end.time) -> tloss[new.ted]
      }
    }
    if (!I[i]) break;
    last.time <- event.times[i]
  }

  return(data.table(
    day = c(event.times[1:i], end.time),
    S = S[c(1:i, i)], I = I[c(1:i, i)], R = R[c(1:i, i)],
    cumulativeI = cumulativeI[c(1:i, i)]
  ))
}

base.het.plot <- function(end.time, maxpop, dt=data.table(runid=integer(), day = numeric(), state=factor(), value=numeric(), runstate=factor())) ggplot() +
  aes(x=day, alpha=runstate, color=state, y=value, group=interaction(runid, state)) +
  scale_alpha_manual(values=c(past=0.3, current=1)) +
  scale_x_continuous(limits = c(0, end.time)) +
  scale_y_continuous(limits = c(0, maxpop)) +
  theme_minimal() + theme(
    axis.title = element_text(size=rel(2)),
    axis.text = element_text(size=rel(2))
  ) +
  guides(alpha="none") + theme(legend.position = "bottom", legend.direction = "horizontal") +
  geom_line(data=dt)


#' @title Sample [het_simulate()]
#'
#' @description
#' Run [het_simulate()] many times with the same parameter combinations.
#'
#' @param n a positive integer, the number of samples to simulate
#'
#' @inheritParams het_simulate
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

  ts <- seq_len(n) |> parallel::mclapply(
    FUN = \(sid, mxdst, params) het_simulate(mxdst, params, sid),
    mxdst = mxdst, params = params,
    mc.preschedule = FALSE
  ) |> rbindlist(idcol = "sample_id")

  return(ts |> setkey(sample_id, day))

}

het.ui <- shinyUI({
  button34 <- actionButton("part34click", "Run")
  fluidPage(
  titlePanel('Heterogeneity Tutorial'),
  navlistPanel(
    tabPanel('Overview', includeMarkdown("inst/hetTut/overview.md"), br()),
    tabPanel('Part 1: Low Variance',
      fluidRow(column(12, includeMarkdown("inst/hetTut/low.md"))),
      fluidRow(
        column(3, numericInput("lowsamples", label = "Add runs:", value = 1, min = 1, max = 50, step = 1)),
        column(3, actionButton("lowclick", "Run")),
        column(3, textOutput("lowruns", inline = T)),
        column(3, textOutput("lowTODO", inline = T))
      ),
      fluidRow(
        column(4, plotOutput("lowhist")), column(4, plotOutput("lowseries")), column(4, plotOutput("lowsizes"))
      ),
      fluidRow(column(12, includeMarkdown("inst/hetTut/more.md"))),
      br()
    ),
    tabPanel('Part 2: Changing Variance and Population',
      fluidRow(column(12, includeMarkdown("inst/hetTut/variance.md"))),
      fluidRow(column(3, numericInput(
          "part3var",
          label = "Variance?",
          value=0.001, min = 0.0001, max = 1
        )),
        column(3, numericInput(
          "part4var",
          label = "Population?",
          value=100, min = 10, max = 500, step = 1
        )),
        column(2, button34),
        column(2, textOutput("part34runs", inline = T)),
        column(2, textOutput("part34TODO", inline = T))
      ),
      fluidRow(
        column(4, plotOutput("part34hist")),
        column(4, plotOutput("part34series")),
        column(4, plotOutput("part34sizes"))
      ),
      fluidRow(

      ),
      br()
    ),
    tabPanel('Part 3: Heterogeneity & R0',
      fluidRow(column(12, includeMarkdown("inst/hetTut/R0.md"))),
      fluidRow(
        column(3,numericInput(
          "part5bmn",
          label = "beta.mean?",
          value=2, min = .1, max = 10
        )),
        column(3,numericInput(
          "part5bvar",
          label = "beta.var?",
          value=.1, min = 1e-6, max = 10
        )),
        column(3,numericInput(
          "part5pop",
          label = "population?",
          value=100, min = 10, max = 500, step = 1
        )),
        column(3,numericInput(
          "part5gmma",
          label = "gamma?",
          value=1, min = 1e-6, max = 10
        ))
      ),
      fluidRow(
        column(3,numericInput(
          "part5rho",
          label = "rho?",
          value=0, min = 0, max = 10
        )),
        column(3, numericInput("part5samples", label = "runs?", value = 10, min = 1, max = 50, step = 1)),
        column(2, actionButton("part5click", "Run")),
        column(2, textOutput("part5runs", inline = T)),
        column(2, textOutput("part5TODO", inline = T))
      ),
      fluidRow(
        column(4, plotOutput("part5hist")),
        column(4, plotOutput("part5series")),
        column(4, plotOutput("part5sizes"))
      ),
      br()
    ), widths = c(2, 10)
  ))
})

emptyseries <- data.table(
  runid=integer(),
  day=numeric(),
  runstate=factor(levels=c("past","current")),
  state=factor(levels=c("S","I","R")),
  value=integer()
)
emptydistro <- data.table(runid=integer(), cumulativeI=integer())

lowmxdst <- het_population(n = 100, beta_mean = 2, beta_var = 0.001)

allowedIterations <- 100

part34samples <- 30

updateSeries <- function(old, new) {
  rbind(
    old[, runstate := "past"],
    new[state != "cumulativeI"]
  )
}

het.server <- shinyServer({

  function(input, output, session) {

    cycleTracking <- reactiveValues(
      lowwas = 0, lowcycles = 0,
      part34was = 0, part34cycles = 0,
      part5was = 0, part5cycles = 0
    )

    rvs <- reactiveValues(
      lowseries = emptyseries,
      lowdistro = emptydistro,
      lowindex = 0,
      part34mxdst = het_population(n = 100, beta_mean = 2, beta_var = 0.001),
      part34distro = emptydistro,
      part34series = emptyseries,
      part34index = 0,
      part5distro = emptydistro,
      part5series = emptyseries,
      part5index = 0
    )

    observe({

      if (input$lowclick - isolate(cycleTracking$lowwas)) {
        cycleTracking$lowcycles <- isolate(input$lowsamples)
        cycleTracking$lowwas <- input$lowclick
      }

      # depend on all buttons
      # lock all buttons

      # do part 12 heavy lifting, if there are cycles available & there's work to do
      isolate(if (cycleTracking$lowcycles) {

          addedts <- het.epidemic.runs(lowmxdst, rvs$lowindex, end.time = 10, gmma = 1)

          rvs$lowseries <- updateSeries(rvs$lowseries, addedts)
          rvs$lowdistro <- rbind(
            rvs$lowdistro,
            addedts[state=="cumulativeI", list(cumulativeI=value[.N]), by=runid]
          )
          # if done, rvs$lowTODO <- FALSE
          rvs$lowindex <- rvs$lowindex + 1
          cycleTracking$lowcycles <- cycleTracking$lowcycles - 1
      })

      if (isolate(cycleTracking$lowcycles)) invalidateLater(0, session)
    })

    observe({
      input$part34click
      rvs$part34index <- 0
      rvs$part34mxdst   <- het_population(n = isolate(input$part4var), beta_mean = 2, beta_var = isolate(input$part3var))
      rvs$part34series  <- emptyseries
      rvs$part34distro  <- emptydistro
    })

    if (!isolate(cycleTracking$lowcycles)) observe({

      if (input$part34click - isolate(cycleTracking$part34was)) {
        cycleTracking$part34cycles <- part34samples
        cycleTracking$part34was <- input$part34click
      }

      isolate(if(cycleTracking$part34cycles) { # do part 34 heavy lifting, if there are cycles available & there's work to do

        addedts <- het.epidemic.runs(rvs$part34mxdst, cycleTracking$part34cycles, end.time = 10, gmma = 1)

        rvs$part34series <- updateSeries(rvs$part34series, addedts)

        rvs$part34distro <- rbind(
          rvs$part34distro,
          addedts[state=="cumulativeI", list(cumulativeI=value[.N]), by=runid]
        )

        rvs$part34index <- rvs$part34index + 1
        cycleTracking$part34cycles <- cycleTracking$part34cycles - 1
      })

      if (cycleTracking$part34cycles) invalidateLater(0, session)

    })

    observe({
      input$part5click
      rvs$part5index <- 0
      rvs$part5mxdst   <- het_population(
        n = isolate(input$part5pop),
        beta_mean = isolate(input$part5bmn),
        beta_var = isolate(input$part5bvar)
      )
      rvs$part5series  <- emptyseries
      rvs$part5distro  <- emptydistro
    })

    if (!isolate(cycleTracking$lowcycles | cycleTracking$part34cycles)) observe({

      if (isolate(cycleTracking$part5was) != input$part5click) {
        cycleTracking$part5cycles <- isolate(input$part5samples)
        cycleTracking$part5was <- input$part5click
      }

      isolate(if(cycleTracking$part5cycles) { # do part 5 heavy lifting, if there are cycles available & there's work to do

        addedts <- het.epidemic.runs(rvs$part5mxdst, cycleTracking$part5cycles, end.time = 10, gmma = input$part5gmma, rho = input$part5rho)
        rvs$part5series <- updateSeries(rvs$part5series, addedts)

        rvs$part5distro <- rbind(
          rvs$part5distro,
          addedts[state=="cumulativeI", list(cumulativeI=value[.N]), by=runid]
        )

        rvs$part5index <- rvs$part5index + 1
        cycleTracking$part5cycles <- cycleTracking$part5cycles - 1
      })

      if (cycleTracking$part5cycles) invalidateLater(0, session)

    })

    output$lowTODO <- renderText(ifelse(cycleTracking$lowcycles, "Running", ""))
    output$lowruns <- renderText(sprintf("Total Runs: %d", rvs$lowindex))

    output$lowhist <- renderPlot(het_plot_hist(lowmxdst, beta_mean = 2))
    output$lowseries <- renderPlot(base.het.plot(10, 100, dt=rvs$lowseries))
    output$lowsizes <- renderPlot({
      if(rvs$lowdistro[,.N != 0]) het.runs.hist(rvs$lowdistro)
    })

    output$part34TODO <- renderText(ifelse(cycleTracking$part34cycles, "Running", ""))
    output$part34runs <- renderText(sprintf("Total Runs: %d", rvs$part34index))

    output$part34hist <- renderPlot(het_plot_hist(rvs$part34mxdst, beta_mean = 2))
    output$part34series <- renderPlot(base.het.plot(10, isolate(input$part4var), dt=rvs$part34series))
    output$part34sizes <- renderPlot({
      if(rvs$part34distro[,.N != 0]) het.runs.hist(rvs$part34distro)
    })

    output$part5TODO <- renderText(ifelse(cycleTracking$part5cycles, "Running", ""))
    output$part5runs <- renderText(sprintf("Total Runs: %d", rvs$part5index))

    output$part5hist <- renderPlot(het_plot_hist(rvs$part5mxdst, beta_mean = isolate(input$part5bmn)))
    output$part5series <- renderPlot(base.het.plot(10, isolate(input$part5pop), dt=rvs$part5series))
    output$part5sizes <- renderPlot({
      if(rvs$part5distro[,.N != 0]) het.runs.hist(rvs$part5distro)
    })

  }})

#' @export
heterogeneityTutorial <- function() shinyApp(ui = het.ui, server = het.server)
