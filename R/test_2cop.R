#' Compute Similarity Between Two Copulas
#'
#' This function performs a statistical test to evaluate the similarity between
#' two copulas represented by the datasets `x` and `y`. It supports both paired
#' and unpaired samples and utilizes Monte Carlo simulation to estimate the
#' p-value and critical value of the test statistic.
#'
#' @param x Numeric matrix representing the first sample of copula data.
#' @param y Numeric matrix representing the second sample of copula data.
#' @param Nsim Integer specifying the number of simulations to perform (default is 100).
#' @param paired Logical indicating whether the samples are paired (default is FALSE).
#' @param alpha Numeric value indicating the confidence level for the test (default is 0.95).
#' 
#' @return A list containing:
#'   - `pvalue`: the p-value of the test indicating similarity.
#'   - `cvm`: the observed test statistic.
#'   - `VaR`: the critical value at the specified confidence level.
#'   - `cvmsim`: the simulated test statistics.
#'
#' @useDynLib TwoCop
#' @importFrom stats quantile
#'
#' @examples
#' x = matrix(runif(100), ncol=2)
#' y = matrix(runif(100), ncol=2)
#' result <- test_2cop(x, y, Nsim=1000, paired=FALSE, alpha=0.95)
#' print(result)
#'
#' @export
test_2cop <- function(x, y, Nsim=100, paired=FALSE, alpha=0.95) {

  dimx  = dim(x)
  dimy =  dim(y)

  n1 = dimx[1]
  n2 = dimy[1]
  d  = dimx[2]
  d2 = dimy[2]

  if (d != d2) stop("Samples x and y must have the same number of dimensions.")
  if (n1 != n2 & paired == TRUE) {
    stop("Paired samples must have the same sample size.")
  }

  out0 = .C("pvalue2",
                  as.double(x),
                  as.double(y),
                  as.integer(n1),
                  as.integer(n2),
                  as.integer(d),
                  as.integer(Nsim),
                  as.integer(paired),
                  cvm    = double(1),
                  cvmhat = double(Nsim),
                  pvalue = double(1),
		  PACKAGE="TwoCop"
                  )

  VaR   = quantile(out0$cvmhat, alpha)
  list(pvalue = out0$pvalue, cvm = out0$cvm,VaR = VaR, cvmsim = out0$cvmhat)
}
