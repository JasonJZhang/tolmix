\name{tolmix-package}
\alias{tolmix-package}
\alias{tolmix}
\docType{package}
\title{Tolerance Intervals for Random and Mixed Effects Models}
\description{
  Tolerance intervals describe a range that contains at least a proportion
  \eqn{P} of a population with confidence \eqn{1-\alpha}. This package
  computes such intervals (and beta-expectation analogues) for the marginal
  distribution of a future observation under variance-component models:

  \itemize{
    \item \code{\link{tolint_oneway}} for the one-way random effects model
      (balanced or unbalanced), with the Mee-Owen (1983) one-sided limits,
      the Hoffman-Kringle (2005) modified large sample two-sided interval,
      and the Krishnamoorthy-Mathew (2004) generalized pivotal quantity
      method.
    \item \code{\link{tolint_lmer}} for general linear mixed models fitted
      with \pkg{lme4}, using a Satterthwaite-type method in the spirit of
      Francq, Lin and Hoyer (2019), a parametric bootstrap pivotal method,
      and beta-expectation (prediction-type) intervals.
  }
}
\keyword{package}
