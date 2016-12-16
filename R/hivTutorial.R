#' @import deSolve data.table ggplot2
NULL

sim.start.year <- 1975
sim.end.year <- 2025
measures <- c("prevalence","incidence","mortality")
measure.cols <- c(prevalence="green", incidence="darkblue", mortality="red")
model.ltys <- c(observed=0, 'model 1'=2, 'model 2'=3, 'model 3'=4)

hiv.pickdataset <- function(ss) hivdata[
  site==ss,
  list(year, value=prevalence, measure=factor("prevalence", levels=measures), name="observed")
]

model.plot <- function(site, sim.start.year, sim.end.year) ggplot() +
  aes(x=year, y=value, color=measure, linetype=name) +
  geom_point(data=hiv.pickdataset(site), size=5, alpha=0.5) +
  theme_minimal() +
  lims(y=c(0,0.50), x=c(sim.start.year, sim.end.year)) +
  scale_color_manual(values=measure.cols, drop=F) +
  scale_linetype_manual(values=model.ltys, drop=T)

invoke.model <- function(
  model, model.xfm, inits, vals, start, end,
  res = 0.1
) melt.data.table(
  model.xfm(data.table(
    lsoda(
      y = inits, times = seq(start, end, res),
      func = model, parms = vals
    )
  )[, year:= time ], vals),
  id.vars="year", measure.vars = measures, variable.name = "measure"
)

model.lines <- function(
  lines=data.table(
    year=numeric(),
    measure=factor(levels=measures),
    value=numeric()
  ), name = "nullmodel"
) geom_line(data=lines[, name := name ])

emptylines <- model.lines()

models <- list(
  `model 1`=list(step=NULL,xfm=NULL)
)

model1xfm <- function(dt, parms) with(parms,
  dt[, prevalence := I/(I+S) ][,
    `:=`(
      incidence = lambda*prevalence*S,
      mortality = delta*prevalence
    )
  ]
)

model1 <- function(t,y,parms) with(c(as.list(y),parms),{
  N <- sum(y)

  inc <- lambda*S*I/N
  dSdt <- birthRate*N - inc - mu*S
  dIdt <- inc - delta*I - mu*I

  list(c(dSdt, dIdt))
})


model2 <- function(t,y,parms) with(c(as.list(y),parms),{
  N <- sum(y)

  lambdaHat <- lambda2*exp(-a*I/N)
  inc <- S*lambdaHat*I/N
  dSdt <- birthRate2*N - inc - mu*S
  dIdt <- inc - delta*I - mu*I

  list(c(dSdt, dIdt))
})

model2xfm <- function(dt, parms) with(parms, {
  dt[, prevalence := I/(I+S) ]
  dt[, incidence := lambda2*exp(-a*prevalence)*prevalence*S ]
  dt[, mortality := delta*prevalence ]
  dt
})

model3 <- function(t,y,parms) with(c(as.list(y),parms),{
  N <- sum(y)
  I <- I1 + I2 + I3 + I4

  lambdaHat <- lambda3*exp(-a3*I/N)
  inc <- S*lambdaHat*I/N

  g <- 4*delta

  dSdt <- birthRate*N - inc - mu*S
  dI1dt <- inc  - (g+mu)*I1
  dI2dt <- g*I1 - (g+mu)*I2
  dI3dt <- g*I2 - (g+mu)*I3
  dI4dt <- g*I3 - (g+mu)*I4

  list(c(dSdt, dI1dt, dI2dt, dI3dt, dI4dt))
})

model3xfm <- function(dt, parms) {
  dt[, I := I1+I2+I3+I4 ]
  dt[, prevalence := I/(I+S) ]
  dt[, incidence := parms$lambda3*exp(-parms$a3*prevalence)*prevalence*S ]
  dt[, mortality := 4*parms$delta*I4/(I+S) ]
  dt
}

make.lines <- function(mk, nm) if (mk) model.lines(invoke.model(

), nm) else emptylines

hiv.server <- shinyServer(
  function(input, output, session){
    refplot <- reactive(model.plot(input$useSite, sim.start.year, sim.end.year))
    model1lines <- reactive(if ("1" %in% input$showmodels) {
      model.lines(invoke.model(
        model1, model1xfm, c(S=1,I=exp(input$p0)), reactiveValuesToList(input),
        sim.start.year, sim.end.year
      ), "model 1")
    } else emptylines)
    model2lines <- reactive(if ("2" %in% input$showmodels) {
      model.lines(invoke.model(
        model2, model2xfm, c(S=1,I=exp(-7)), modifyList(reactiveValuesToList(input), list(mu = 0.018, delta = 0.1)),
        sim.start.year, sim.end.year
      ), "model 2")
    } else emptylines)
    model3lines <- reactive(if ("3" %in% input$showmodels) {
      model.lines(invoke.model(
        model3, model3xfm, c(S=1,I1=exp(-7), I2=0,I3=0,I4=0),
        modifyList(reactiveValuesToList(input), list(birthRate = 0.029, mu = 0.018, delta = 0.1)),
        sim.start.year, sim.end.year
      ), "model 3")
    } else emptylines)
    output$modelplot <- renderPlot(refplot()+model1lines()+model2lines()+model3lines())
  }
)

#' @export
hivTutorial <- function() shinyApp(ui = shinyUI(
  withMathJax(fluidPage(
    sidebarLayout(
      sidebarPanel(
        h3('HIV in Harare tutorial'),
        withMathJax(includeMarkdown("inst/hiv/sidebar.md")),
        selectInput("useSite", "Dataset for fitting:",
          choices = hivdata[, levels(site)],
          selected = 'Harare', width = '200px'
        ),
        checkboxGroupInput("showmodels", "Show Model(s):", choices=1:3, inline = T)
      ),
      mainPanel(
        plotOutput("modelplot", height="300px", width = "600px"),
        tabsetPanel(
        tabPanel("Model 1", includeMarkdown("inst/hiv/model1.md"),
          # actionButton("runButton1","Show model output (see plot above)", width = "100%"),
          sliderInput('lambda',
            div(HTML('Rate at which new infections occur (&lambda;):')),
            min = 0, step = 0.05, max = 0.6,
            value = 0.5
          ),
          sliderInput('p0',
            'Natural log of the intial infection prevalence:',
            min = -10, step = 0.5, max = -3,
            value = -7
          ),
          sliderInput('birthRate',
            'Per capita birth rate (b):',
            min = 0.01, step = 0.001, max = 0.04,
            value = 0.029
          ),
          sliderInput('delta',
            div(HTML('Per capita disease-induced mortality rate (&delta;):')),
            min = 0, step = 0.05, max = 0.5,
            value = 0.1
          ),
          sliderInput('mu',
            div(HTML('Per capita background mortality rate (&mu;):')),
            min = 0.01, step = 0.001, max = 0.04,
            value = 0.018
          )
        ),
        tabPanel("Model 2", includeMarkdown("inst/hiv/model2.md"),
          # actionButton("runButton2","Show model output (see plot above)", width = "100%"),
          sliderInput('a',
            'Rate of decline as a function of prevalence (a):',
            min = 0, step = 0.5, max = 20,
            value = 0
          ),
          sliderInput('lambda2',
            div(HTML('Rate at which new infections occur (&lambda;):')),
            min = 0, step = 0.05, max = 0.6,
            value = 0.5
          ),
          sliderInput('birthRate2',
            div(HTML('Per capita birth rate (b):')),
            min = 0.01, step = 0.001, max = 0.04,
            value = 0.029
          )
        ),
        tabPanel("Model 3", includeMarkdown("inst/hiv/model3.md"),
          # actionButton("runButton3","Show model output (see plot above)",width = "100%"),
          sliderInput('a3',
            div(HTML('Rate of decline as a function of prevalence (a):')),
            min = 0, step = 0.5, max = 20,
            value = 8
          ),
          sliderInput('lambda3',
            div(HTML('Rate at which new infections occur (&lambda;):')),
            min = 0, step = 0.05, max = 0.6,
            value = 0.5
          )
        )
      ), style = "height:600px; overflow-y:scroll;")
    )
  ))
), server = hiv.server)
