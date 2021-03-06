% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tail.R
\name{tailNormal}
\alias{tailNormal}
\title{Normal approximation for tail estimation}
\usage{
tailNormal(x, dfg1, dfg2)
}
\arguments{
\item{x}{points to evaluate tail probabilities in}

\item{dfg1}{dfg object specifying null model}

\item{dfg2}{dfg object specifying foreground model}
}
\value{
A dataframe with columns, x, tail estimate and confidence intervals
}
\description{
Normal approximation for tail estimation
}

