% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM tmm.R
\name{generate_modt}
\alias{generate_modt}
\title{Generate moderated t-statistics from the neural network t-mixture model}
\usage{
generate_modt(fit)
}
\arguments{
\item{fit}{list returned by \link[modtcov]{fit.t.model}.}
}
\description{
Generate moderated t-statistics from the neural network t-mixture model
}
\examples{
mod.t.cov.df <- create.t.stat.cov.df(exampleY, exampleX, c(1,-1,0,0,0))
fit <- fit.t.model(mod.t.cov.df)
generate_modt(fit)
}
