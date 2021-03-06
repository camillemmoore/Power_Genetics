% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood_calc_function.R
\name{ll.linear.selector}
\alias{ll.linear.selector}
\title{Function to return log likelihood function for specified model type}
\usage{
ll.linear.selector(model)
}
\arguments{
\item{model}{The genetic model in the linear regression: "Dominant", "Additive", "Recessive", "2df" or "null"}
}
\value{
Log likelihood function for specified model type
}
\description{
Convenience function to return log likelihood function for specified model type
}
\examples{
ll.linear.selector("Dominant")

}
