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

het.run <- function(rid, mxdst, end.time, gmma, rho) {
  set.seed(rid)
  pop.size <- length(mxdst)
  tloss <- trec <- tinf <- rep(NaN, pop.size)

  new.ted <- sample.int(pop.size, 1, prob = mxdst)
  new.time <- 0 # first individual is infected at time 0
  (( new.time -> tinf[new.ted] ) +
      rexp(1, gmma) -> trec[new.ted] ) +
        suppressWarnings(rexp(1, rho)) -> tloss[new.ted]
  ## initialize time series data frame
  ts <- data.table(day=new.time, S=pop.size-1, I=1, R = 0, cumulativeI=1)
  next.time <- data.table::copy(ts)
  ## calculate current inf rate - should this be changing?  was inside while loop
  inf.rate <- mean(mxdst) * pop.size

  ## expected samples to get sum(samples) > end.time
  # times <- rexp(ceiling(end.time*1.1*inf.rate), inf.rate)
  # while (sum(times) < end.time) times <- c(times, rexp(ceiling(end.time*0.1*inf.rate), inf.rate))

  while(next.time[, (day < end.time) & (I > 0)]) {
    ## Choose potential event time
    new.time <- rexp(1, inf.rate) + next.time$day
    ## figure out if anyone has recovered between last event and this one; I->R
    recoveries <- sum((trec > next.time$day) & (trec < new.time), na.rm = T)
    next.time[,`:=`(I=I-recoveries, R=R+recoveries)]
    ## figure out if anyone has lost immunity between last event and this one; R->S
    who.loss <- which((tloss > next.time$day) & (tloss < new.time))
    ilosses <- length(who.loss)
    if (ilosses > 0) { # if there are immunity losses...
      ## reset tinf, trec, tloss for individuals who lost immunity
      tinf[who.loss] <- trec[who.loss] <- tinf[who.loss] <- NaN
      ## update R, S counts
      next.time[,`:=`(R=R-ilosses, S=S+ilosses)]
    }

    ## choose infec*tor*
    cur.tor <- sample(1:pop.size, 1, prob = mxdst)
    ## has this id been infected?
    cur.inf <- {
      tr <- trec[cur.tor]
      !is.nan(tr) & (tr > new.time)
    }
    # if source infected, then choose (potentially) infec*ted* individual
    if(cur.inf) {
      # choose from everyone except infec*tor*
      new.ted <- sample((1:pop.size)[-cur.tor], 1, prob = mxdst[-cur.tor])
      if(is.nan(tinf[new.ted])) { ## is this ID currently susceptible?
        ## then infect them
        next.time[,`:=`(
          cumulativeI = cumulativeI + 1,
          S = S - 1, I = I + 1
        )]
        # infected @ new.time, recover @ infected + recovery draw, loss @ recover + loss draw
        (( new.time -> tinf[new.ted] ) +
            rexp(1, gmma) -> trec[new.ted] ) +
              suppressWarnings(rexp(1, rho)) -> tloss[new.ted]
      }
    }
    ts <- rbind(ts, next.time)
    next.time <- data.table::copy(next.time)[, day := new.time ]
  }
  return(ts[, runid := rid ])
}

base.het.plot <- function(end.time, maxpop, dt=data.table(runid=integer(), day = numeric(), state=factor(), value=numeric())) ggplot() +
  aes(x=day, alpha=runid, color=state, y=value, group=interaction(runid, state)) +
  scale_alpha_continuous(limits = c(-10, NA)) +
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
  ts <- het.run(runs[1], mxdst, end.time, gmma, rho)
  if (!is.null(tell.run)) tell.run$inc(1/length(runs), detail = paste("Finished run", runs[1]))
  for (runid in runs[-1]) {
      ts <- rbind(ts, het.run(runid, mxdst, end.time, gmma, rho))
      if (!is.null(tell.run)) tell.run$inc(1/length(runs), detail = paste("Finished run", runid))
  }
  return(melt.data.table(ts, id.vars = c("runid","day"), variable.name = "state"))
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
    tabPanel('Parts 1 & 2: Low Variance',
      fluidRow(column(12, includeMarkdown("inst/hetTut/part1.md"))),
      fluidRow(
        column(3, numericInput("part1samples", label = "Add runs:", value = 1, min = 1, max = 50, step = 1)),
        column(3, actionButton("part1click", "Run")),
        column(3, textOutput("part1runs", inline = T)),
        column(3, textOutput("part1TODO", inline = T))
      ),
      fluidRow(
        column(4, plotOutput("part1hist")), column(4, plotOutput("part1series")), column(4, plotOutput("part1sizes"))
      ),
      fluidRow(column(12, includeMarkdown("inst/hetTut/part2.md"))),
      br()
    ),
    tabPanel('Part 3: Changing Variance and Population',
      fluidRow(column(12, includeMarkdown("inst/hetTut/part3.md"))),
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
    tabPanel('Part 4: Heterogeneity & R0',
      fluidRow(column(12, includeMarkdown("inst/hetTut/part4.md"))),
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
  state=factor(levels=c("S","I","R")),
  value=integer()
)
emptydistro <- data.table(runid=integer(), cumulativeI=integer())

part1mxdst <- het.population(pop.size = 100, beta.mean = 2, beta.var = 0.001)

allowedIterations <- 100

part34samples <- 30

het.server <- shinyServer({

  function(input, output, session) {

    cycleTracking <- reactiveValues(
      part1was = 0, part1cycles = 0,
      part34was = 0, part34cycles = 0,
      part5was = 0, part5cycles = 0
    )

    rvs <- reactiveValues(
      part1series = emptyseries,
      part1distro = emptydistro,
      part1index = 0,
      part34mxdst = het.population(pop.size = 100, beta.mean = 2, beta.var = 0.001),
      part34distro = emptydistro,
      part34series = emptyseries,
      part34index = 0,
      part5distro = emptydistro,
      part5series = emptyseries,
      part5index = 0
    )

    observe({

      if (input$part1click - isolate(cycleTracking$part1was)) {
        cycleTracking$part1cycles <- isolate(input$part1samples)
        cycleTracking$part1was <- input$part1click
      }

      # depend on all buttons
      # lock all buttons

      # do part 12 heavy lifting, if there are cycles available & there's work to do
      isolate(if (cycleTracking$part1cycles) {

          addedts <- het.epidemic.runs(part1mxdst, rvs$part1index, end.time = 10, gmma = 1)

          rvs$part1series <- rbind(
            rvs$part1series,
            addedts[state != "cumulativeI"]
          )
          rvs$part1distro <- rbind(
            rvs$part1distro,
            addedts[state=="cumulativeI", list(cumulativeI=value[.N]), by=runid]
          )
          # if done, rvs$part1TODO <- FALSE
          rvs$part1index <- rvs$part1index + 1
          cycleTracking$part1cycles <- cycleTracking$part1cycles - 1
      })

      if (isolate(cycleTracking$part1cycles)) invalidateLater(0, session)
    })

    observe({
      input$part34click
      rvs$part34index <- 0
      rvs$part34mxdst   <- het.population(pop.size = isolate(input$part4var), beta.mean = 2, beta.var = isolate(input$part3var))
      rvs$part34series  <- emptyseries
      rvs$part34distro  <- emptydistro
    })

    if (!isolate(cycleTracking$part1cycles)) observe({

      if (input$part34click - isolate(cycleTracking$part34was)) {
        cycleTracking$part34cycles <- part34samples
        cycleTracking$part34was <- input$part34click
      }

      isolate(if(cycleTracking$part34cycles) { # do part 34 heavy lifting, if there are cycles available & there's work to do

        addedts <- het.epidemic.runs(rvs$part34mxdst, cycleTracking$part34cycles, end.time = 10, gmma = 1)

        rvs$part34series <- rbind(
          rvs$part34series,
          addedts[state != "cumulativeI"]
        )

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

    if (!isolate(cycleTracking$part1cycles | cycleTracking$part34cycles)) observe({

      if (isolate(cycleTracking$part5was) != input$part5click) {
        cycleTracking$part5cycles <- isolate(input$part5samples)
        cycleTracking$part5was <- input$part5click
      }

      isolate(if(cycleTracking$part5cycles) { # do part 5 heavy lifting, if there are cycles available & there's work to do

        addedts <- het.epidemic.runs(rvs$part5mxdst, cycleTracking$part5cycles, end.time = 10, gmma = input$part5gmma, rho = input$part5rho)

        rvs$part5series <- rbind(
          rvs$part5series,
          addedts[state != "cumulativeI"]
        )

        rvs$part5distro <- rbind(
          rvs$part5distro,
          addedts[state=="cumulativeI", list(cumulativeI=value[.N]), by=runid]
        )

        rvs$part5index <- rvs$part5index + 1
        cycleTracking$part5cycles <- cycleTracking$part5cycles - 1
      })

      if (cycleTracking$part5cycles) invalidateLater(0, session)

    })

    output$part1TODO <- renderText(ifelse(cycleTracking$part1cycles, "Running", ""))
    output$part1runs <- renderText(sprintf("Total Runs: %d", rvs$part1index))

    output$part1hist <- renderPlot(het.hist(part1mxdst, beta.mean = 2))
    output$part1series <- renderPlot(base.het.plot(10, 100, dt=rvs$part1series))
    output$part1sizes <- renderPlot({
      if(rvs$part1distro[,.N != 0]) het.runs.hist(rvs$part1distro)
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
