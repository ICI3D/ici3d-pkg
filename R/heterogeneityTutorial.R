#' @import data.table ggplot2 shiny
NULL

het.population <- function(
  pop.size, beta.mean, beta.var, seed = 0
) {
  set.seed(seed)
  theta = beta.var / beta.mean
  kk = beta.mean / theta
  mxdst <- rgamma(pop.size, shape = kk, scale = theta)
  mxdst[rev(order(mxdst))] ## order individuals by riskiness
}

het.hist <- function(mxdst, beta.mean, wd = .1) {
  breaks <- seq(0, ceiling(max(mxdst)), by = wd)
  p <- ggplot(data.frame(het = mxdst)) +
    geom_histogram(
      aes(x=het, y=..density..),
      breaks = breaks
    ) + xlab("contact rate (1/day)") +
    geom_vline(xintercept=beta.mean, color="red") +
    theme_bw() + theme(
      panel.border = element_blank(),
      axis.title = element_text(size=rel(2)),
      axis.text = element_text(size=rel(2))
    )
  p
}

# het.run <- function(rid, mxdst, end.time, gmma, rho) {
#   set.seed(rid)
#   pop.size <- length(mxdst)
#   tloss <- trec <- tinf <- rep(NaN, pop.size)
#
#   new.ted <- sample.int(pop.size, 1, prob = mxdst)
#   new.time <- 0 # first individual is infected at time 0
#   (( new.time -> tinf[new.ted] ) +
#       rexp(1, gmma) -> trec[new.ted] ) +
#         suppressWarnings(rexp(1, rho)) -> tloss[new.ted]
#   ## initialize time series data frame
#   ts <- data.table(day=new.time, S=pop.size-1, I=1, R = 0, cumulativeI=1)
#   next.time <- data.table::copy(ts)
#   ## calculate current inf rate - should this be changing?  was inside while loop
#   inf.rate <- mean(mxdst) * pop.size
#
#   ## expected samples to get sum(samples) > end.time
#   # times <- rexp(ceiling(end.time*1.1*inf.rate), inf.rate)
#   # while (sum(times) < end.time) times <- c(times, rexp(ceiling(end.time*0.1*inf.rate), inf.rate))
#
#   while(next.time[, (day < end.time) & (I > 0)]) {
#     ## Choose potential event time
#     new.time <- rexp(1, inf.rate) + next.time$day
#     ## figure out if anyone has recovered between last event and this one; I->R
#     recoveries <- sum((trec > next.time$day) & (trec < new.time), na.rm = T)
#     next.time[,`:=`(I=I-recoveries, R=R+recoveries)]
#     ## figure out if anyone has lost immunity between last event and this one; R->S
#     who.loss <- which((tloss > next.time$day) & (tloss < new.time))
#     ilosses <- length(who.loss)
#     if (ilosses > 0) { # if there are immunity losses...
#       ## reset tinf, trec, tloss for individuals who lost immunity
#       tinf[who.loss] <- trec[who.loss] <- tinf[who.loss] <- NaN
#       ## update R, S counts
#       next.time[,`:=`(R=R-ilosses, S=S+ilosses)]
#     }
#
#     ## choose infec*tor*
#     cur.tor <- sample(1:pop.size, 1, prob = mxdst)
#     ## has this id been infected?
#     cur.inf <- {
#       tr <- trec[cur.tor]
#       !is.nan(tr) & (tr > new.time)
#     }
#     # if source infected, then choose (potentially) infec*ted* individual
#     if(cur.inf) {
#       # choose from everyone except infec*tor*
#       new.ted <- sample((1:pop.size)[-cur.tor], 1, prob = mxdst[-cur.tor])
#       if(is.nan(tinf[new.ted])) { ## is this ID currently susceptible?
#         ## then infect them
#         next.time[,`:=`(
#           cumulativeI = cumulativeI + 1,
#           S = S - 1, I = I + 1
#         )]
#         # infected @ new.time, recover @ infected + recovery draw, loss @ recover + loss draw
#         (( new.time -> tinf[new.ted] ) +
#             rexp(1, gmma) -> trec[new.ted] ) +
#               suppressWarnings(rexp(1, rho)) -> tloss[new.ted]
#       }
#     }
#     ts <- rbind(ts, next.time)
#     next.time <- data.table::copy(next.time)[, day := new.time ]
#   }
#   return(ts[, runid := rid ])
# }

het.run <- function(rid, mxdst, end.time, gmma, rho) {
  set.seed(rid)
  pop.size <- length(mxdst)
  tloss <- trec <- rep(NaN, pop.size)

  inf.rate <- mean(mxdst) * pop.size

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

  return(data.table(runid = rid, day = event.times[1:i], S=S[1:i], I=I[1:i], R=R[1:i], cumulativeI = cumulativeI[1:i]))
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
  mxdst = het.population(pop.size = 300, beta.mean = 1, beta.var = .5),
  runs,
  gmma = 1, rho = 0, # lose immunity
  end.time = 5, tell.run = NULL
) {
  ts <- het.run(runs[1], mxdst, end.time, gmma, rho)[, runstate := "current"]
  if (!is.null(tell.run)) tell.run$inc(1/length(runs), detail = paste("Finished run", runs[1]))
  for (runid in runs[-1]) {
      ts <- rbind(ts[, runstate := "past"], het.run(runid, mxdst, end.time, gmma, rho)[, runstate := "current"])
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

lowmxdst <- het.population(pop.size = 100, beta.mean = 2, beta.var = 0.001)

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
      part34mxdst = het.population(pop.size = 100, beta.mean = 2, beta.var = 0.001),
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
      rvs$part34mxdst   <- het.population(pop.size = isolate(input$part4var), beta.mean = 2, beta.var = isolate(input$part3var))
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
        pop.size = isolate(input$part5pop),
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
