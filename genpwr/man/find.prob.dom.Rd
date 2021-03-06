% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odds_ratio_function.R
\name{find.prob.dom}
\alias{find.prob.dom}
\title{Dominant probability finding function}
\usage{
find.prob.dom(x, P_AA, P_AB, P_BB, cr, like)
}
\arguments{
\item{x}{Probability of disease given AB or BB}

\item{P_AA}{Probability the allele is homozygous for the major allele}

\item{P_AB}{Probability the allele is heterozygous}

\item{P_BB}{Probability the allele is homozygous for the minor allele}

\item{cr}{proportion of cases in the sample (cases/(cases + controls)).}

\item{like}{Expected log likelihood}
}
\value{
: The "a" in the binomial function ax^2 + bx + c that arises in solution for the additive OR functions
}
\description{
Operates within add.or.function to find probability of disease in a dominant truth given AB or BB, additive test model
}
\examples{
find.prob.dom(0.1510677, 0.5625, 0.375, 0.0625, 0.3, -0.57162)

}
