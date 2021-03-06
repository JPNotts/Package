% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GaussianProcess.R
\name{get_l_single}
\alias{get_l_single}
\title{Approximate l from PWC}
\usage{
get_l_single(
  x.pwc,
  end.time.real,
  end.time = 20,
  sigma.f2 = 500,
  sigma.n2 = 1e-05
)
}
\arguments{
\item{x.pwc}{the intensity function from a PWC prior.}

\item{end.time.real}{The end.time of the experiment}

\item{end.time}{End time used in the statistical analysis.}

\item{sigma.f2}{sigma.f2}

\item{sigma.n2}{sigma.n2}
}
\value{
l value.
}
\description{
Approximate l from PWC
}
