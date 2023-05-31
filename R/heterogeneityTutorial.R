#' @import data.table ggplot2 shiny
NULL

#' @title Sample Population with Heterogeneous Contact Rates
#'
#' @description Create a sample of individuals with heterogeneous contact rates
#'
#' @param n Number of individuals in the population
#'
#' @param beta.mean Mean contact rate
#'
#' @param beta.var Variance of contact rate
#'
#' @param seed Seed for random number generation
#'
#' @return A numeric vector of contact rates, sorted from highest to lowest
#'
#' @details Per `?rgamma`, the gamma distribution is parameterized by `shape`
#' and `scale` parameters.  These parameters are related to the mean by
#' \deqn{E(X)=shape*scale} and the variance by \deqn{Var(X)=shape*scale^2}.
#'
#' @family heterogeneity
#' @export
#' @examples
#' require(ICI3D)
#' egpop <- het.population(1000, 2, .5)
#' het.hist(egpop, 2)
het.population <- function(
  n, beta.mean, beta.var, seed
) {
  if (!missing(seed)) { set.seed(seed) }
  scale <- beta.var / beta.mean
  shape <- beta.mean / scale
  sort(rgamma(n, shape = shape, scale = scale), decreasing = TRUE)
}

#' @title Plot Histogram of Contact Rates
#'
#' @description Plot a histogram of contact rates, with annotations for
#' population and sample means.
#'
#' @param mxdst A numeric vector of contact rates, as produced by
#' [het.population()].
#'
#' @param beta.mean The mean contact rate, as passed to [het.population()].
#'
#' @inheritParams ggplot2::geom_histogram
#' @inheritDotParams ggplot2::geom_histogram
#'
#' @return A `ggplot` object
#' 
#' @family heterogeneity
#' @export
het.hist <- function(
  mxdst, beta.mean, binwidth = .1, ...
) {
  return(ggplot(data.table(het = mxdst)) +
    aes(x = het) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = binwidth, ...
    ) +
    geom_vline(xintercept = beta.mean, color = "red") +
    geom_blank(data = function(dt) { dt[, {
      r <- range(het)
      .(het = c(floor(r[1]), ceiling(r[2])))
    }] }) +
    labs(x = "contact rate (1/day)") +
    theme_bw() + theme(
      panel.border = element_blank(),
      axis.title = element_text(size = rel(2)),
      axis.text = element_text(size = rel(2))
  ))
}

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
  return(list(source, sink))
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
#' @return A numeric vector of event times: the cumulative sum of draws from
#' `sampler` until the last value exceeds `maxval`. That sum is then trimmed to
#' values `<= maxval`. The first value in the vector is always 0.
#'
#' @family heterogeneity
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


#' @title Generate an SIR Heterogeneous Mixing Sample
#'
#' @param seed Random number seed
#'
#' @param mxdst A numeric vector of contact rates, as produced by
#' [het.population()].
#'
#' @param tmax The maximum time to simulate until
#'
#' @param gmma The mean recovery rate (I -> R transition rate)
#'
#' @param rho The mean immunity loss rate (R -> S transition rate)
#'
#' @return A data.table
#' 
#' @family heterogeneity
#' @export 
#' @examples
#' het.run(
#'   mxdst = het.population(n = 100, beta.mean = 2, beta.var = 0.001),
#'   tmax = 10, gmma = 1
#' )
het.run <- function(
  mxdst, tmax, gmma, rho = 0, seed = NULL
) {
  # if using a seed ...
  if (!is.null(seed)) { set.seed(seed) }

  # get the population size
  pop.size <- length(mxdst)
  # initialize vectors for individual loss / recovery times
  # NaN => happens at time infinity, i.e. never
  tloss <- trec <- rep(NaN, pop.size)

  # calculate the average time-to-interaction
  contact.rate <- sum(mxdst)

  # pre-draw all the interaction times
  event.times <- sample_until(
    tmax, 1/contact.rate, rexp, rate = contact.rate
  )

  # define draw functions for recovery and immnunity-loss delays
  rec_delays <- function() { rexp(1, gmma) }
  # if rho is 0, then we don't need to draw loss delays
  loss_delays <- if (rho > 0) {
    function() { rexp(1, rho) }
  } else {
    function() { tmax }
  }

  # pre-allocate outcome vectors
  stepmax <- length(event.times)
  curI <- 1 # number of currently infected - if ever < 1, epidemic is over
  # delX: change in X at this time
  delS <- c(pop.size - curI, rep.int(0, stepmax - 1))
  delI <- c(curI, rep.int(0, stepmax - 1))
  delR <- rep.int(0, stepmax)
  # cumX: cumulative X at this time
  cumI <- delI

  contacts <- setNames(sample_pair(stepmax, mxdst), c("contactor", "contactee"))

  # who was initially infected? when do they recover / lose immunity?
  init.tee <- contacts$contactor[1]
  last.time <- event.times[1]
  trec[init.tee] <- last.time + rec_delays()
  tloss[init.tee] <- trec[init.tee] + loss_delays()

  # for all other events ...
  for (i in 2:stepmax) {
    # update the time
    last.time <- event.times[i]
    # check for recoveries since last event
    who.recover <- which(trec < last.time)
    if (recoveries <- length(who.recover)) {
      trec[who.recover] <- NaN
      curI <- curI - recoveries
      delI[i] <- -recoveries
      delR[i] <- recoveries
    }

    # check for immunity losses since last event
    who.loss <- which(tloss < last.time)
    if (ilosses <- length(who.loss)) {
      tloss[who.loss] <- NaN
      delR[i] <- -ilosses
      delS[i] <- ilosses
    }

    cur.tor <- contacts$contactor[i]
    cur.tee <- contacts$contactee[i]
     # if X has a recovery time == currently infected
     # and Y has no recovery / loss time == currently susceptible

    ## is this individual infected?
    if (!(is.nan(tloss[cur.tor]) || is.nan(trec[cur.tor]))) {
      if (is.nan(tloss[cur.tee])) { ## is the contactee currently susceptible?
        ## then infect them
        delS[i] <- delS[i] - 1
        delI[i] <- delI[i] + 1
        cumI[i] <- 1

        trec[cur.tee] <- last.time + rec_delays()
        tloss[cur.tee] <- trec[cur.tee] + loss_delays()

        curI <- curI + 1
      }
    }
    if (curI == 0) break;
  }

  return(data.table(
    runid = seed,
    day = event.times[1:i],
    S = cumsum(delS[1:i]),
    I = cumsum(delI[1:i]),
    R = cumsum(delR[1:i]),
    cumulativeI = cumsum(cumI[1:i])
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


het.epidemic.runs <- function(
  mxdst = het.population(n = 300, beta.mean = 1, beta.var = .5),
  runs,
  gmma = 1, rho = 0, # lose immunity
  end.time = 5, tell.run = NULL
) {
  ts <- het.run(mxdst, end.time, gmma, rho, runs[1])[, runstate := "current"]
  if (!is.null(tell.run)) tell.run$inc(1/length(runs), detail = paste("Finished run", runs[1]))
  for (runid in runs[-1]) {
      ts <- rbind(ts[, runstate := "past"], het.run(mxdst, end.time, gmma, rho, runid)[, runstate := "current"])
      if (!is.null(tell.run)) tell.run$inc(1/length(runs), detail = paste("Finished run", runid))
  }
  return(melt.data.table(ts, id.vars = c("runid","day","runstate"), variable.name = "state"))
}

het.runs.hist <- function(ts) ggplot(
  ts[,list(`final size` = max(cumulativeI)), by=runid]
) + aes(x=`final size`) + geom_histogram(binwidth = 5) +
  theme_minimal() + theme(
    axis.title = element_text(size=rel(2)),
    axis.text = element_text(size=rel(2))
  ) + ylab("# of outbreaks")

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

lowmxdst <- het.population(n = 100, beta.mean = 2, beta.var = 0.001)

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
      part34mxdst = het.population(n = 100, beta.mean = 2, beta.var = 0.001),
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
      rvs$part34mxdst   <- het.population(n = isolate(input$part4var), beta.mean = 2, beta.var = isolate(input$part3var))
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
      rvs$part5mxdst   <- het.population(
        n = isolate(input$part5pop),
        beta.mean = isolate(input$part5bmn),
        beta.var = isolate(input$part5bvar)
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

    output$lowhist <- renderPlot(het.hist(lowmxdst, beta.mean = 2))
    output$lowseries <- renderPlot(base.het.plot(10, 100, dt=rvs$lowseries))
    output$lowsizes <- renderPlot({
      if(rvs$lowdistro[,.N != 0]) het.runs.hist(rvs$lowdistro)
    })

    output$part34TODO <- renderText(ifelse(cycleTracking$part34cycles, "Running", ""))
    output$part34runs <- renderText(sprintf("Total Runs: %d", rvs$part34index))

    output$part34hist <- renderPlot(het.hist(rvs$part34mxdst, beta.mean = 2))
    output$part34series <- renderPlot(base.het.plot(10, isolate(input$part4var), dt=rvs$part34series))
    output$part34sizes <- renderPlot({
      if(rvs$part34distro[,.N != 0]) het.runs.hist(rvs$part34distro)
    })

    output$part5TODO <- renderText(ifelse(cycleTracking$part5cycles, "Running", ""))
    output$part5runs <- renderText(sprintf("Total Runs: %d", rvs$part5index))

    output$part5hist <- renderPlot(het.hist(rvs$part5mxdst, beta.mean = isolate(input$part5bmn)))
    output$part5series <- renderPlot(base.het.plot(10, isolate(input$part5pop), dt=rvs$part5series))
    output$part5sizes <- renderPlot({
      if(rvs$part5distro[,.N != 0]) het.runs.hist(rvs$part5distro)
    })

  }})

#' @export
heterogeneityTutorial <- function() shinyApp(ui = het.ui, server = het.server)
