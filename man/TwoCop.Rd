\name{TwoCop}
\alias{TwoCop}

\title{Nonparametric test of equality between two copulas
}
\description{This function performs the nonparametric test of equality
between two copulas proposed by Remillard and Scaillet (2009). The test is based on the Cramer-von-Mises statistic between the two
empirical copulas. An approximate p-value is returned.
}
\usage{
TwoCop(x, y, Nsim=100, paired=FALSE, alpha=0.95)
}

\arguments{
  \item{x}{
 \code{n} by \code{d} matrix containing the first dataset.
}
  \item{y}{ 
 \code{m} by \code{d} matrix containing the second dataset.
}
  \item{Nsim}{
  Number of iterations used in the approximation of the p-value.

}
  \item{paired}{\code{FALSE} (default) means that x and y are from two independent
  populations, \code{TRUE} indicates paired data.
}
  \item{alpha}{Level of the calculated VaR. Default is 0.95.
}
}
\details{
Details of the method can be found in Remillard and Scaillet (2009).
}
\value{ A list of the following objects:
  \item{cvm}{Value of the Cramer-von Mises test statistic.}
  \item{pvalue}{pvalue based on the multiplier Monte Carlo method 
         with \code{Nsim} iterations.}
  \item{cvmsim}{Simulated values of the Cramer-von Mises statistic.}
  \item{VaR}{\code{alpha} quantile of the simulated Cramer-von Mises statistics.}
}
\references{
Remillard, B. & Scaillet, O. (2009) Testing for equality between two copulas. Journal of Multivariate Analysis, 100, 377-386.
}
\author{Bruno Remillard and Jean-Francois Plante

}

\examples{

# Simulating a bivariate normal (copula = independence)

  X=matrix(rnorm(100),ncol=2)

# Simulating a bivriate exponential distribution with a Clayton copula

  v=runif(50)
  theta=1
  x<-1/(1/runif(50)/v^(theta+1))^(1/(theta+1))
  u<-(x^(-theta)-v^(-theta)+1)^(-1/theta)
  Y=cbind(-log(1-u),-log(1-v))
  
# Testing equality of the copulas

  TwoCop(X,Y)$pvalue
}
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ htest }
\keyword{ copula }
\keyword{ ranks }
