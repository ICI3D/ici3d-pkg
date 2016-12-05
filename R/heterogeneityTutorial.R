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
  p <- ggplot(data.frame(heterogeniety = mxdst)) +
    geom_histogram(
      aes(x=heterogeniety, y=..density..),
      breaks = breaks
    ) +
    geom_vline(xintercept=beta.mean, color="red") +
    theme_bw() + theme(panel.border = element_blank())
  p
}

het.run <- function(runid, mxdst, end.time, gmma, rho) {
  set.seed(runid)
  pop.size <- length(mxdst)
  tloss <- trec <- tinf <- rep(NaN, pop.size)

  new.ted <- sample.int(pop.size, 1, prob = mxdst)
  new.time <- 0 # first individual is infected at time 0
  (( new.time -> tinf[new.ted] ) +
      rexp(1, gmma) -> trec[new.ted] ) +
        rexp(1, rho) -> tloss[new.ted]
  ## initialize time series data frame
  ts <- data.table(time=new.time, S=pop.size-1, I=1, R = 0, cumulativeI=1)
  next.time <- data.table::copy(ts)
  ## calculate current inf rate - should this be changing?  was inside while loop
  inf.rate <- mean(mxdst) * pop.size

  ## expected samples to get sum(samples) > end.time
  # times <- rexp(ceiling(end.time*1.1*inf.rate), inf.rate)
  # while (sum(times) < end.time) times <- c(times, rexp(ceiling(end.time*0.1*inf.rate), inf.rate))

  while(next.time[, (time < end.time) & (I > 0)]) {
    ## Choose potential event time
    new.time <- rexp(1, inf.rate) + next.time$time
    ## figure out if anyone has recovered between last event and this one; I->R
    recoveries <- sum((trec > next.time$time) & (trec < new.time), na.rm = T)
    next.time[,`:=`(I=I-recoveries, R=R+recoveries)]
    ## figure out if anyone has lost immunity between last event and this one; R->S
    who.loss <- which((tloss > next.time$time) & (tloss < new.time))
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
              rexp(1, rho) -> tloss[new.ted]
      }
    }
    ts <- rbind(ts, next.time)
    next.time <- data.table::copy(next.time)[, time := new.time ]
  }
  return(ts[, runid := runid ])
}

base.het.plot <- function(end.time, maxpop, dt=data.table(runid=integer(), time = numeric(), state=factor(), value=numeric())) ggplot() +
  aes(x=time, alpha=runid, color=state, y=value, group=interaction(runid, state)) +
  scale_alpha_continuous(limits = c(-10, NA)) +
  scale_x_continuous(limits = c(0, end.time)) +
  scale_y_continuous(limits = c(0, maxpop)) +
  theme_minimal() + guides(alpha="none") + theme(legend.position = "bottom", legend.direction = "horizontal") +
  geom_line(data=dt)


het.epidemic.runs <- function(
  # runs = 1,
  mxdst = het.population(pop.size = 300, beta.mean = 1, beta.var = .5),
  runs,
  gmma = 1, rho = 1/10, # lose immunity
  end.time = 5, tell.run = NULL
) {
  ts <- het.run(runs[1], mxdst, end.time, gmma, rho)
  if (!is.null(tell.run)) tell.run$inc(1/length(runs), detail = paste("Finished run", runs[1]))
  for (runid in runs[-1]) {
      ts <- rbind(ts, het.run(runid, mxdst, end.time, gmma, rho))
      if (!is.null(tell.run)) tell.run$inc(1/length(runs), detail = paste("Finished run", runid))
  }
  return(melt.data.table(ts, id.vars = c("runid","time"), variable.name = "state"))
}

het.runs.hist <- function(ts) ggplot(
  ts[,list(fsize = max(cumulativeI)), by=runid]
) + aes(x=fsize) + geom_histogram(binwidth = 5) +
  theme_minimal() + xlim(0, NA)
# hist(ts[,list(fsize = max(cumulativeI)),by=runid]$fsize, breaks = 100, xlab = "cumulative # infected (final size)",
#      xlim = c(0, f.size),
#      ylab = "frequency", main = "outbreak size distribution", col = "black")


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
    tabPanel('Part 4: Heterogeniety & R0',
      fluidRow(column(12, includeMarkdown("inst/hetTut/part4.md"))),
      fluidRow(
        column(3,numericInput(
          "part5var",
          label = "R0?",
          value=2, min = .1, max = 10
        )),
        column(3, actionButton("part5click", "Run")),
        column(3, textOutput("part5runs", inline = T)),
        column(3, textOutput("part5TODO", inline = T))
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


# want to update graphics between work
# server block needs to a maximum amount of work (aka, simulation runs)
# then release for render updates, but resume after if there's still work to do

emptyseries <- data.table(runid=integer(), time=numeric(), state=factor(levels=c("S","I","R")), value=integer())
emptydistro <- data.table(runid=integer(), cumulativeI=integer())


het.server <- shinyServer({

  part1mxdst <- het.population(pop.size = 100, beta.mean = 2, beta.var = 0.001)

  allowedIterations <- 100

  part34samples <- 30

  cycleTracking <- reactiveValues(
    part1was = 0, part1cycles = 0,
    part34was = 0, part34cycles = 0,
    part5was = 0, part5TODO = F, part5cycles = 0,
    itersLeft = allowedIterations
  )

  rvs <- reactiveValues(
    part1series = emptyseries,
    part1distro = emptydistro,
    part1index = 0,
    part34mxdst = het.population(pop.size = 100, beta.mean = 2, beta.var = 0.001),
    part34distro = emptydistro,
    part34series = emptyseries,
    part34index = 0
  )

  function(input, output, session) {

    observe({
      input$part5click
      rvs$part5Amxdst <- het.population(
        pop.size = 100, beta.mean = isolate(input$part5var), beta.var = 0.1
      )
      rvs$part5Bmxdst <- het.population(
        pop.size = 100, beta.mean = isolate(input$part5var), beta.var = 10
      )
      rvs$part5Aseries <- base.het.plot(end.time = 10, maxpop=100)
      rvs$part5Bseries <- base.het.plot(end.time = 10, maxpop=100)
    })

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

    # observe({
    #
    #   if (isolate(cycleTracking$part5was) != input$part5click) {
    #     cycleTracking$part5TODO <- TRUE
    #     cycleTracking$part5cycles <- part34samples
    #     cycleTracking$part5was <- input$part5click
    #   }
    #
    #   isolate(if(cycleTracking$part5TODO & (cycleTracking$itersLeft > 0)) { # do part 5 heavy lifting, if there are cycles available & there's work to do
    #
    #     cycleTracking$part5cycles <- cycleTracking$part5cycles - 1
    #     cycleTracking$part5TODO <- cycleTracking$part5cycles != 0
    #     cycleTracking$itersLeft <- cycleTracking$itersLeft - 1
    #   })
    #
    # })

#    cycleTracking$itersLeft <- allowedIterations

#    print(sprintf("updating iters left: %d", isolate(cycleTracking$itersLeft)))
    #    output$part1cycles <- renderText(isolate(cycleTracking$part1cycles))
    #    output$itersLeft <- renderText(isolate(cycleTracking$itersLeft))


    output$part1TODO <- renderText(ifelse(cycleTracking$part1cycles, "Running", "Waiting"))
    output$part1runs <- renderText(sprintf("Total Runs: %d", rvs$part1index))

    output$part1hist <- renderPlot(het.hist(part1mxdst, beta.mean = 2))
    output$part1series <- renderPlot(base.het.plot(10, 100, dt=rvs$part1series))
    output$part1sizes <- renderPlot({
      if(rvs$part1distro[,.N != 0]) het.runs.hist(rvs$part1distro)
    })

    output$part34TODO <- renderText(ifelse(cycleTracking$part34cycles, "Running", "Waiting"))
    output$part34runs <- renderText(sprintf("Total Runs: %d", rvs$part34index))

    output$part34hist <- renderPlot(het.hist(rvs$part34mxdst, beta.mean = 2))
    output$part34series <- renderPlot(base.het.plot(10, isolate(input$part4var), dt=rvs$part34series))
    output$part34sizes <- renderPlot({
      if(rvs$part34distro[,.N != 0]) het.runs.hist(rvs$part34distro)
    })

    # output$part5Ahist <- renderPlot(het.hist(rvs$part5Amxdst, beta.mean = input$part5var))
    # output$part5Aseries <- renderPlot({ rvs$part5Aseries })
    # output$part5Asizes <- renderPlot({
    #   if(rvs$part5Adistro[,.N != 0]) het.runs.hist(rvs$part5Adistro)
    # })
    # output$part5Bhist <- renderPlot(het.hist(rvs$part5Bmxdst, beta.mean = input$part5var))
    # output$part5Bseries <- renderPlot({ rvs$part5Bseries })
    # output$part5Bsizes <- renderPlot({
    #   if(rvs$part5Bdistro[,.N != 0]) het.runs.hist(rvs$part5Bdistro)
    # })

  }})

#' @export
heterogeneityTutorial <- function() shinyApp(ui = het.ui, server = het.server)
