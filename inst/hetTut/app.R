#' @import data.table ggplot2 shiny
NULL

het.population <- function(
  pop.size, beta.mean, beta.var
) {
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

base.het.plot <- function(end.time, maxpop) ggplot() +
  aes(x=time, alpha=factor(runid), color=state, y=value) +
  scale_alpha_discrete(limits = c(-10, NA)) +
  scale_x_continuous(limits = c(0, end.time)) +
  scale_y_continuous(limits = c(0, maxpop)) +
  theme_minimal() + guides(alpha="none") +
  geom_line(data=data.table(runid=integer(), time = numeric(), state=factor(), value=numeric()))

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
) + aes(x=fsize) + geom_histogram()
# hist(ts[,list(fsize = max(cumulativeI)),by=runid]$fsize, breaks = 100, xlab = "cumulative # infected (final size)",
#      xlim = c(0, f.size),
#      ylab = "frequency", main = "outbreak size distribution", col = "black")


ui <- shinyUI({
fluidPage(
  titlePanel('Dynamical Fever: computer exercise'),
  navlistPanel(
    tabPanel('Overview', includeMarkdown("overview.md")),
    tabPanel('Parts 1 & 2: Low Variance',
      includeMarkdown("part1.md"),
      numericInput("part1samples", label = "Run how many additional simulations?", value = 1, min = 1, max = 50, step = 1), actionButton("part1click", "Run Simulations"),
      plotOutput("part1hist"),
      plotOutput("part1series"),
      plotOutput("part1sizes"),
      includeMarkdown("part2.md")
    ),
    tabPanel('Parts 3 & 4: Changing Variance and Population',
      includeMarkdown("part3.md"),
      numericInput(
        "part3var",
        label = "Desired variance?",
        value=0.001, min = 0.0001, max = 1
      ),
      actionButton("part34click", "Run Simulations"),
      plotOutput("part34hist"),
      plotOutput("part34series"),
      includeMarkdown("part4.md"),
      numericInput(
        "part4var",
        label = "Population size?",
        value=100, min = 10, max = 500, step = 1
      ),
      actionButton("part34click", "Run Simulations")
    ),
    tabPanel('Part 5: Heterogeniety & R0',
      includeMarkdown("part5.md"),
      actionButton("part5click", "Run Simulations"),
      numericInput(
        "part5var",
        label = "R0?",
        value=2, min = .1, max = 10
      ),
      plotOutput("part5hist"),
      plotOutput("part5series")
    )
  ))
})


server <- shiny::shinyServer({
  part1mxdst <- het.population(pop.size = 100, beta.mean = 2, beta.var = 0.001)
  rvs <- reactiveValues(
    part1series = base.het.plot(end.time = 10, maxpop=100),
    part1distro = data.table(runid=integer(), cumulativeI=integer()),
    part1index = 0
  )
  function(input, output) {
    output$part1hist <- renderPlot(het.hist(part1mxdst, beta.mean = 2))

    # part 1
    observe({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      clicks <- input$part1click
      freezeReactiveValue(input, "part1click")
      start <- isolate(rvs$part1index) + 1
      end <- start - 1 + isolate(input$part1samples)
      rvs$part1index <- end
      progress$inc(0, detail=sprintf("starting %d samples...", isolate(input$part1samples)))
      addedts <- het.epidemic.runs(part1mxdst, start:end, end.time = 10, gmma = 1, tell.run = progress)

      rvs$part1series <- isolate(rvs$part1series) + geom_line(data=addedts[state != "cumulativeI"])
      rvs$part1distro <- rbind(
        isolate(rvs$part1distro),
        addedts[state=="cumulativeI", list(cumulativeI=value[.N]), by=runid]
      )
    })

    # part 2
    # observe({
    #   runid <- input$part2click
    #   freezeReactiveValue(input, "part2click")
    #   rvs$part1series <- isolate(rvs$part1series) + het.epidemic.runs(part1mxdst, runid, end.time = 10, gmma = 1)
    #   rvs$part1samps <- runid
    # })

    output$part1series <- renderPlot({ rvs$part1series })
    output$part1sizes <- renderPlot({
      if(rvs$part1distro[,.N != 0]) het.runs.hist(rvs$part1distro)
    })

  }})

shinyApp(ui = ui, server = server)

# het.epidemic(beta.mean = 2, beta.var = .001, runs = 5, end.time = 10, pop.size = 100, gmma = 1)
