% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulatingSpikes.R
\name{simulate_spikes}
\alias{simulate_spikes}
\title{Simulate spike sequences}
\usage{
simulate_spikes(
  end.time,
  int.fn,
  hyper,
  steps = 1000,
  T.min = NULL,
  ISI.type = "Gamma",
  sequences = 1,
  add.end = TRUE,
  do.log = T
)
}
\arguments{
\item{end.time}{The length of time we want to simulate spikes}

\item{int.fn}{The intensity function}

\item{hyper}{The ISI parameter}

\item{steps}{The number of steps to split the experiment into}

\item{T.min}{Refractory period (Where pdf = 0)}

\item{ISI.type}{The ISI distribution}

\item{sequences}{The number of spike sequences to generate}

\item{add.end}{Do we include the end time on the end of the spikes}

\item{do.log}{Flag for whether we do the calculations on the log scale}
}
\value{
Matrix containing multiple spikes sequences
}
\description{
Simulate spike sequences
}
