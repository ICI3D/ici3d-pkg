#' @import data.table ggplot2 shiny
NULL

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

rf.cols <- c("black","red","blue","green","yellow","purple")

rf.ui <- app.ui.generic('Reed-Frost Exploration',
  tabPanel('Overview', includeMarkdown("inst/rf/overview.md"), br()),
  tabPanel('Deterministic Reed-Frost',
    includeMarkdown("inst/rf/det.md"),
    numericInput("detN", "Population", value=100, min=10, max=10000, step=1),
    sliderInput("detp", "Infection Probability", value=0.1, min = 0, max = 0.25, step = 0.001),
    plotOutput("detpanel"),
    textOutput("detfsize"),
    br()
  ),
  tabPanel('Stochastic Reed-Frost', includeMarkdown("inst/rf/stoch.md"),
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
))


# want to update graphics between work
# server block needs to a maximum amount of work (aka, simulation runs)
# then release for render updates, but resume after if there's still work to do

emptystoch <- data.table(
  runid = integer(),
  time = integer(),
  color = factor(levels=rf.cols),
  variable = factor(levels=c("S","C","R")),
  count = integer()
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
          id.vars = c("runid", "time", "color"), value.name = "count"
        )
      )
    })

    output$detpanel <- renderPlot({
      ggplot(melt.data.table(rvs$detrun, id.vars = "time", value.name = "count")) +
        theme_bw() + theme(strip.text = element_text(size = rel(5))) +
        aes(x=time, y=count) + facet_grid(. ~ variable) + geom_step()
    })
    output$detfsize <- renderText(sprintf("final size: %f", rvs$detrun[,sum(C)]))

    output$stochpanel <- renderPlot({
      ggplot(rvs$stochruns) + theme_bw() + theme(strip.text = element_text(size = rel(5))) +
        aes(x=time, y=count, group=runid, color=factor(color, levels=rf.cols)) +
        facet_grid(. ~ variable) + geom_step() +
        scale_color_manual(values = rf.cols, drop=F)
    })

    output$stochfsize <- renderPlot({
      distro <- rvs$stochruns[variable=="C",list(`total cases`=sum(count)),by=list(runid, color)]
      ggplot(distro) + theme_minimal() +
        aes(x=`total cases`, fill=factor(color, levels=rf.cols)) + geom_bar() +
        scale_fill_manual(values = rf.cols, drop=F)
    })

}})

#' @export
reedFrost <- function() shinyApp(ui = rf.ui, server = rf.server)
