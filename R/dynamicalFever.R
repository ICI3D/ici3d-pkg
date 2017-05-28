#' @import data.table ggplot2 shiny

epi.duration <- function(epi=run.example()) with(epi, {
  c(
    Dogs=diff(range(Time[Cases.Pop1>0]))+1,
    People=diff(range(Time[Cases.Pop2>0]))+1
  )
})

plot.cases <- function(epi=run.example()){
  barplot(epi$Cases.Pop1,names.arg=epi$Time,
          cex.lab=2,
          ylim=c(0,200),
          cex.axis = 1.1,
          xlab='Time',
          ylab='Cases',
          cex.main = 1.5,
          main='Dogs',
          xaxt='n') -> bb
  axis(1,bb[seq(0,30,5)+1],seq(0,30,5),cex.axis = 1.1,xpd=T,line=.2)
  barplot(epi$Cases.Pop2,names.arg=epi$Time,
          cex.lab=2,
          ylim=c(0,200),
          cex.axis = 1.1,
          xlab='Time',
          ylab='Cases',
          cex.main = 1.5,
          main='People',
          xaxt='n') -> bb
  axis(1,bb[seq(0,30,5)+1],seq(0,30,5),cex.axis = 1.1,xpd=T,line=.2)
}

plot.example <- function(epi=run.example(),plot.Re=FALSE){
  if(plot.Re){
    par(mar=c(5,5,2,1),mfrow=c(3,2))
  }else{
    par(mar=c(5,5,2,1),mfrow=c(2,2))
  }
  plot.cases(epi)
  plot(epi$Time,epi$FOI.Pop1,
       type='s',      # Use a 'step' plot because time is treated as discrete
       bty='n',
       lwd=2,
       cex.lab=2,
       ylim=c(0,1),
       cex.axis = 1.1,
       xlab='Time',
       ylab='Force of Infection',
       cex.main = 1.5,
       main='')
  plot(epi$Time,epi$FOI.Pop2,
       type='s',      # Use a 'step' plot because time is treated as discrete
       bty='n',
       lwd=2,
       cex.lab=2,
       ylim=c(0,1),
       cex.axis = 1.1,
       xlab='Time',
       ylab='Force of Infection',
       cex.main = 1.5,
       main='')
  if(plot.Re){
    plot(epi$Time,epi$FOI.Pop1*epi$Sus.Pop1/epi$Cases.Pop1,
         type='p',      # Use a 'step' plot because time is treated as discrete
         bty='n',
         lwd=2,
         col=2,
         cex.lab=2,
         ylim=c(0,2),
         cex.axis = 1.1,
         xlab='Time',
         ylab=expression(R[effective]),
         cex.main = 1.5,
         main='')
    abline(h=1)
    lines(epi$Time,2/999*epi$Sus.Pop1,col=4)
  }
}

run.example <- function(
  VaxPct.Pop1=0, VaxPct.Pop2=0
) sir.ex1.cb(VaxPct.Pop1, VaxPct.Pop2)

make.pop <- function(N, vaxRate, I) {
  stopifnot(is.integer(N), is.integer(I), 0 < N, I <= N, 0 <= vaxRate, vaxRate <= 1)
  as.list(environment())
}

sim.dynFever <- function(
  dog.pop = make.pop(1000, 0, 1), human.pop = make.pop(1000, 0, 0),
  R0=2, MAXTIME=30, pp=0.002
) {
  qq.dog <- 1-R0/(dog.pop$N-1)
  qq.human <- 1-pp
  res <- data.table(
    time=1:MAXTIME,
    cases.dog=0,   sus.dog=0,   foi.dog=0,
    cases.human=0, sus.human=0, foi.human=0
  )
  cases.dog <- integer(MAXTIME)
  cases.dog[1] <- dog.pop$N
  foi.dog <- numeric(MAXTIME)

  return(res)
}

sir.ex1.cb <- function(
  VaxPct.Pop1=0, VaxPct.Pop2=0,
  R0=2, N1=1000, N2=1000,
  MAXTIME=30, I1.0=1,
  pp=0.002
){

  qq.1 <- 1-R0/(N1-1)  # Pairwise probability of avoiding potentially infectious contact (Pop1)
  qq.2 <- 1-pp        # Pairwise probability of avoiding potentially infectious contact (Pop2)

  Cases.Pop1 <- I1.0
  Sus.Pop1 <- max(0,N1-I1.0-round(VaxPct.Pop1*N1/100))
  FOI.Pop1 <- (1-qq.1^Cases.Pop1)

  Cases.Pop2 <- 0
  Sus.Pop2 <- N2-round(VaxPct.Pop2*N2/100)
  FOI.Pop2 <- (1-qq.2^Cases.Pop1)

  for(Time in 1:MAXTIME){

    Cases.Pop1 <- c(Cases.Pop1,rbinom(1,Sus.Pop1[Time],FOI.Pop1[Time]))
    Sus.Pop1 <- c(Sus.Pop1,Sus.Pop1[Time]-Cases.Pop1[Time+1])
    FOI.Pop1 <- c(FOI.Pop1,(1-qq.1^Cases.Pop1[Time+1]))

    Cases.Pop2 <- c(Cases.Pop2,rbinom(1,Sus.Pop2[Time],FOI.Pop2[Time]))
    Sus.Pop2 <- c(Sus.Pop2,Sus.Pop2[Time]-Cases.Pop2[Time+1])
    FOI.Pop2 <- c(FOI.Pop2,(1-qq.2^Cases.Pop1[Time+1]))
  }
  return(data.frame(
    Time=0:MAXTIME,
    Cases.Pop1,
    Sus.Pop1,
    Cases.Pop2,
    Sus.Pop2,
    FOI.Pop1,
    FOI.Pop2
  ))
}

total.cases <- function(epi=run.example()){
  c(Dogs=sum(epi$Cases.Pop1),People=sum(epi$Cases.Pop2))
}

vax.eff <- function(VAXPCT,POP,REPS=100){
  switch(POP,
         Pop1 = {
           replicate(REPS,sum(run.example(VaxPct.Pop1=VAXPCT)$Cases.Pop1))
         },
         Pop2 = {
           replicate(REPS,sum(run.example(VaxPct.Pop2=VAXPCT)$Cases.Pop2))
         }
  )
}

df.ui <- shinyUI(fluidPage(
  titlePanel('Dynamical Fever: computer exercise'),
  navlistPanel(
    tabPanel('Overview', includeMarkdown("inst/dynFev/overview.md"), br()),
    tabPanel('Part 1: Epidemic dynamics', includeMarkdown("inst/dynFev/part1.md"), br()),
    tabPanel('Part 2: Introduction of a veterinary vaccine', includeMarkdown("inst/dynFev/part2.md"), br()),
    tabPanel('Part 3: Introduction of a human vaccine', includeMarkdown("inst/dynFev/part3.md"), br()),
    tabPanel('Part 4: Moving forward',
      h1('Part 4: Moving forward'),
      p('Decide on target levels of vaccination for dogs and people in 2017, keeping in mind that it is unlikely that you will be able to acheive 100 percent vaccination of either population. Enter these values below, each as a number between 0 and 100.'),
      sliderInput('VaxPct.Dogs',
        'Target vaccination level for DOG population:',
        min = 0, max = 100, value = 0,
        width='350px'
      ),
      sliderInput('VaxPct.Humans',
        'Target vaccination level for HUMAN population:',
        min = 0, max = 100, value = 0,
        width='350px'
      ),
      #),
      p("We'll now run the model once to see an example of what might happen if these levels of vaccination were acheived in 2017. Scroll down to see what happened..."),
      plotOutput('targetPlot'),
      p(strong('Is this what you expected to happen? You can reload the page as many times as you like to get a feeling for whether the outcome above is typical of what would be expected when these levels of vaccination are achieved.')),
      p('Now let\'s run the simulation 1000 times with the target vaccination levels. This may take a while.'),
      p('These results can now be plotted to give you a better feeling for the variation in outcomes under an intervention acheiving the targeted levels of vaccination in each population:'),
      plotOutput('distPlot'),
      p(strong('Do these plots for your chosen target vaccination levels give you any additional insight into the processes underlying DF transmission? If not, try lowering your target vaccination levels for at least one of the populations and repeating this section. What is each of these plots showing, and do the results surprise you?')),
      br()
    ),
    tabPanel('Part 5: Vaccination outcomes', includeMarkdown("inst/dynFev/part5.md"), br())
  ))
)

df.server <- shinyServer({
  function(input, output) {
    output$targetPlot <- renderPlot({
      target.2017 <- run.example(input$VaxPct.Dogs, input$VaxPct.Humans)
      par(mar=c(5,5,5,1), mfrow=c(1,2)) # Set up plot
      plot.cases(target.2017)
    })
    output$distPlot <- renderPlot({
      target.runs <- replicate(1000,total.cases(run.example(input$VaxPct.Dogs,input$VaxPct.Humans)))
      par(mar=c(5,5,5,1),mfrow=c(1,2)) # Set up plot
      hist(target.runs['Dogs',], col='dark grey',
           main='Dogs',
           xlab='Number of canine cases',
           ylab='Number of runs')

      hist(target.runs['People',], col='dark grey',
           main='People',
           xlab='Number of human cases',
           ylab='Number of runs')
      #     plot(target.runs['Dogs',],target.runs['People',],
      #          main='For each of 1000 runs',
      #          xlab='Number of canine cases',
      #          ylab='Number of human cases')
    })

  }})

#' @export
dynamicalFever <- function() shinyApp(ui = df.ui, server = df.server)
