% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eda.R
\name{eda_plot}
\alias{eda_plot}
\title{ggPlot Diagnostic Plots}
\usage{
eda_plot(y, type = c("run", "lag", "hist", "qqn", "spectrum"))
}
\arguments{
\item{y}{Numeric vector.}

\item{type}{Which plot to run.  Default is "run".}
}
\description{
This function takes in a vector of numerical data and a type of ggplot
in order to inspect with more detail.
}
\examples{
eda_ggPlot(rnorm(500))
eda_ggPlot(rnorm(500), type = "hist")
}
