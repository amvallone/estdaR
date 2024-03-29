% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theta.R
\name{theta}
\alias{theta}
\title{Rank Depcompostion index}
\usage{
theta(x, Regime, nsim = NULL)
}
\arguments{
\item{x}{a \eqn{nxk} matrix with \eqn{l \geq 2} successive columns of a variable are later moments in time}

\item{Regime}{values corresponding to which regime each observation belongs to}

\item{nsim}{number of random spatial permutations to generate for computationally based inference}
}
\value{
A data frame
}
\description{
Compute the Regime mobility measure Rey (2004)
}
\details{
For sequence of time periods Theta measures the extent to which rank changes for a variable measured over n locations are in the same direction within mutually exclusive and exhaustive partitions (regimes) of the n locations.
}
\examples{
data(mexico)
theta(mexico[,1:7],mexico[,8],999)

}
\references{
Rey, S.J. (2004) “Spatial dependence in the evolution of regional income distributions,” in A. Getis, J. Mur and H.Zoeller (eds). Spatial Econometrics and Spatial Statistics. Palgrave, London, pp. 194-213.
}
