% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/acceptance_functions.R
\name{boltzmannAccept}
\alias{boltzmannAccept}
\title{Probabilistic acceptance based on Boltzmann's distribution}
\usage{
boltzmannAccept(delta, temperature = 1, ...)
}
\arguments{
\item{delta}{Increment in the evaluation function computed as the evaluation of the new solution minus the evaluation of the existing solution}

\item{temperature}{Temperature parameter in Boltzmann's distribution. This temperature has to be a value strictly greater than 0. By default this parameter is equal to 1}

\item{...}{This parameter is ignored}
}
\value{
\code{TRUE} with probability given by exp(-1*delta/tempreature)
}
\description{
This function accepts solutions (returns \code{TRUE}) with probability given by Boltzmann's distribution
}
\examples{
boltzmannAccept (delta=-1, temperature=2)
boltzmannAccept (delta=30, temperature=200)
boltzmannAccept (delta=30, temperature=2)

}
\seealso{
Other Acceptance functions: \code{\link{alwaysAccept}},
  \code{\link{thresholdAccept}}
}

