% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heterogeneityTutorial.R
\name{het.run}
\alias{het.run}
\title{Generate an SIR Heterogeneous Mixing Sample}
\usage{
het.run(mxdst, tmax, gmma, rho = 0, seed = NULL)
}
\arguments{
\item{mxdst}{A numeric vector of contact rates, as produced by
\code{\link[=het.population]{het.population()}}.}

\item{tmax}{The maximum time to simulate until}

\item{gmma}{The mean recovery rate (I -> R transition rate)}

\item{rho}{The mean immunity loss rate (R -> S transition rate)}

\item{seed}{Random number seed}
}
\value{
A data.table
}
\description{
Generate an SIR Heterogeneous Mixing Sample
}
\examples{
het.run(
  mxdst = het.population(n = 100, beta.mean = 2, beta.var = 0.001),
  tmax = 10, gmma = 1
)
}
\seealso{
Other heterogeneity: 
\code{\link{het.hist}()},
\code{\link{het.population}()},
\code{\link{sample_pair}()},
\code{\link{sample_until}()}
}
\concept{heterogeneity}
