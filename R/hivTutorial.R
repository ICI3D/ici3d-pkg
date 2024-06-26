#' @import deSolve
NULL

#' @export
hivTutorial <- function() rmarkdown::run(system.file("hiv", "index.Rmd", package = "ICI3D"))

#' @title HIV in Harare Models
#'
#' @description
#' These are the five models considered in the [hivtutorial()] interactive
#' shiny app. They all conform to the requirements for [deSolve::ode()]
#' function arguments.
#'
#' @details
#'  - `hiv_SI`: a susceptible-infectious model, with no other processes
#'  - `hiv_SIhet`: `hiv_SI`, but with less transmission as more of the
#'  population is infectious, representing higher-risk people being infected
#'  first
#'  - `hiv_SI4`: `hiv_SIhet`, but with the infectious compartment represented
#'  a boxcar
#'  - `hiv_SI4control`: `hiv_SI4`, but introducing a control program
#'  - `hiv_SI4react`: `hiv_SI4`, but with individuals reacting to AIDS deaths
#'
#' @rdname hararemodels
#' @export
hiv_SI <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    # computed states
    N <- sum(y)

    # processes
    birth <- b * N
    infection <- lambda * S * I / N
    death <- mu * y
    disease_death <- delta * I

    # states
    dSdt <- +(birth) - (infection)
    dIdt <- +(infection) - (disease_death)

    return(list(
      c(dSdt, dIdt) - death # all must die
    ))
  })
}

#' @rdname hararemodels
#' @export
hiv_SIhet <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    # computed states
    N <- sum(y)

    # adjustment: represent depletion of high-risk susceptibles
    lambdaHat <- lambda * exp(-a * I / N)

    # processes
    birth <- b * N
    infection <- lambdaHat * S * I / N
    death <- mu * y
    disease_death <- delta * I

    # states
    dSdt <- +(birth) - (infection)
    dIdt <- +(infection) - (disease_death)

    return(list(
      c(dSdt, dIdt) - death # all must die
    ))
  })
}

#' @rdname hararemodels
#' @export
hiv_SI4 <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    # computed states
    N <- sum(y)
    I <- I1 + I2 + I3 + I4

    # adjustment: represent depletion of high-risk susceptibles
    lambdaHat <- lambda * exp(-a * I / N)
    # adjustment: boxcars need to have the same ultimate leaving rate
    g <- 4 * delta

    # processes
    birth <- b * N
    infection <- lambdaHat * S * I / N
    death <- mu * y
    progression <- g * c(I1, I2, I3, I4)

    # states
    dSdt <- +(birth) - (infection)
    dIdt <- +(c(infection, progression[1:3])) - (progression)

    return(list(
      c(dSdt, dIdt) - death # all must die
    ))
  })
}

#' @rdname hararemodels
#' @export
hiv_SI4control <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    # computed states
    N <- sum(y)
    I <- I1 + I2 + I3 + I4

    # adjustment: represent depletion of high-risk susceptibles
    lambdaHat <- lambda * exp(-a * I / N)
    # adjustment: boxcars need to have the same ultimate leaving rate
    g <- 4 * delta
    # adjustment: control program
    control_effect <- min(1, 1 - cMax / (1 + exp(-cRate * (t - cHalf))))

    # processes
    birth <- b * N
    infection <- control_effect * lambdaHat * S * I / N
    death <- mu * y
    progression <- g * c(I1, I2, I3, I4)

    # states
    dSdt <- +(birth) - (infection)
    dIdt <- +(c(infection, progression[1:3])) - (progression)

    return(list(
      c(dSdt, dIdt) - death # all must die
    ))
  })
}

#' @rdname hararemodels
#' @export
hiv_SI4react <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    # computed states
    N <- sum(y)
    I <- I1 + I2 + I3 + I4

    # adjustment: represent depletion of high-risk susceptibles
    lambdaHat <- lambda * exp(-a * I / N)
    # adjustment: boxcars need to have the same ultimate leaving rate
    g <- 4 * delta
    # adjustment: reaction to AIDS deaths
    death_reaction <- exp(-q * g * I4 / N)

    # processes
    birth <- b * N
    infection <- death_reaction * lambdaHat * S * I / N
    death <- mu * y
    progression <- g * c(I1, I2, I3, I4)

    # states
    dSdt <- +(birth) - (infection)
    dIdt <- +(c(infection, progression[1:3])) - (progression)

    return(list(
      c(dSdt, dIdt) - death # all must die
    ))
  })
}
