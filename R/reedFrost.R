#' @import data.table ggplot2 shiny
NULL

rf.hist <- function(dt) {
  p <- ggplot(dt) +
    geom_histogram(
      aes(x=fsize, y=..density..),
      breaks = breaks
    ) +
    geom_vline(xintercept=beta.mean, color="red") +
    theme_bw() + theme(panel.border = element_blank())
  p
}

det.rf.run <- function(N, p, minC = 0.1) {
  result <- data.table(S=N-1, C=1, R=0)
  nC <- 1
  while(nC > minC) {
    nC <- result[.N, S*(1-(1-p)^C)]
    result <- rbind(
      result,
      result[.N, list(S=S-nC, C=nC, R=R+C)]
    )
  }
  result[, time := .I-1 ]
}

sto.rf.run <- function(N, p) {
  result <- data.table(S=N-1, C=1, R=0)
  nC <- 1
  while(nC) {
    nC <- result[.N, rbinom(1,S,(1-(1-p)^C))]
    result <- rbind(
      result,
      result[.N, list(S=S-nC, C=nC, R=R+C)]
    )
  }
  result[, time := .I-1 ]
}

rf.hist <- function(ts) ggplot(
  ts[,list(fsize = max(cumulativeI)), by=runid]
) + aes(x=fsize) + geom_histogram(binwidth = 5) +
  theme_minimal() + xlim(0, NA)
# hist(ts[,list(fsize = max(cumulativeI)),by=runid]$fsize, breaks = 100, xlab = "cumulative # infected (final size)",
#      xlim = c(0, f.size),
#      ylab = "frequency", main = "outbreak size distribution", col = "black")

rf.cols <- c("black","red","blue","green","yellow","purple")

rf.ui <- shinyUI({
  fluidPage(
    titlePanel('Reed-Frost Exploration'),
    navlistPanel(
      tabPanel('Overview', includeMarkdown("inst/rf/overview.md"), br()),
      tabPanel('Deterministic Reed-Frost',
        includeMarkdown("inst/rf/det.md"),
        numericInput("detN", "Population", value=100, min=10, max=10000, step=1),
        sliderInput("detp", "Infection Probability", value=0.1, min = 0, max = 0.25, step = 0.001),
        plotOutput("detpanel"),
        textOutput("detfsize"),
        br()
      ),
      tabPanel('Stochastic Reed-Frost',
        includeMarkdown("inst/rf/stoch.md"),
        numericInput("stochN", "Population", value=100, min=10, max=10000, step=1),
        sliderInput("stochp", "Infection Probability", value=0.1, min = 0, max = 1, step = 0.001),
        selectInput("stochcol", "line color",
          choices=rf.cols,
          selected = "black"
        ),
        actionButton("runstoch","Run"),
        actionButton("resetstoch","Reset"),
        plotOutput("stochpanel"),
        plotOutput("stochfsize"),
        br()
      ),
      widths = c(2, 10)
    )
  )
})


# want to update graphics between work
# server block needs to a maximum amount of work (aka, simulation runs)
# then release for render updates, but resume after if there's still work to do

emptyseries <- data.table(runid=integer(), time=numeric(), state=factor(levels=c("S","I","R")), value=integer())
emptydistro <- data.table(runid=integer(), cumulativeI=integer())

emptystoch <- data.table(
  runid = integer(),
  time = integer(),
  color = factor(levels=rf.cols),
  variable = factor(levels=c("S","C","R")),
  value = integer()
)

rf.server <- shinyServer({
  rvs <- reactiveValues(
    detrun = det.rf.run(100, 0.1),
    stochruns = emptystoch
  )
  function(input, output, session) {
    observe({
      rvs$detrun <- det.rf.run(input$detN, input$detp)
    })

    observe({
      input$resetstoch
      rvs$stochruns <- emptystoch
    })

    observe({
      input$runstoch
      rid <- suppressWarnings(isolate(rvs$stochruns)[.N, max(max(runid)+1,1)])
      set.seed(rid)
      rvs$stochruns <- rbind(
        isolate(rvs$stochruns),
        melt.data.table(
          sto.rf.run(isolate(input$stochN), isolate(input$stochp))[, runid := rid ][, color := isolate(input$stochcol) ],
          id.vars = c("runid", "time", "color")
        )
      )
    })

    output$detpanel <- renderPlot({
      ggplot(melt.data.table(rvs$detrun, id.vars = "time")) +
        theme_minimal() +
        aes(x=time, y=value) + facet_grid(. ~ variable) + geom_step()
    })
    output$detfsize <- renderText(sprintf("final size: %f", rvs$detrun[,sum(C)]))

    output$stochpanel <- renderPlot({
      ggplot(rvs$stochruns) + theme_minimal() +
        aes(x=time, y=value, group=runid, color=factor(color, levels=rf.cols)) +
        facet_grid(. ~ variable) + geom_step() +
        scale_color_manual(values = rf.cols, drop=F)
    })

    output$stochfsize <- renderPlot({
      distro <- rvs$stochruns[variable=="C",list(fsize=sum(value)),by=list(runid, color)]
      ggplot(distro) + theme_minimal() +
        aes(x=fsize, fill=factor(color, levels=rf.cols)) + geom_bar() +
        scale_fill_manual(values = rf.cols, drop=F)
    })

}})

#' @export
reedFrost <- function() shinyApp(ui = rf.ui, server = rf.server)
