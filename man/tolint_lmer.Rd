\name{tolint_lmer}
\alias{tolint_lmer}
\title{Tolerance intervals for linear mixed models fitted with lme4}
\description{
  Computes one-sided or two-sided tolerance intervals for the distribution
  of a future observation under a linear mixed-effects model fitted with
  \code{lme4::lmer()}. The future observation is taken at the fixed-effect
  covariates given in \code{newdata} and is assumed to come from \emph{new}
  levels of every grouping factor, so its marginal distribution is
  \deqn{Y \sim N(x_0'\beta, \; \sigma_T^2), \qquad
        \sigma_T^2 = \sum_j z_{0j}' G_j z_{0j} + \sigma_e^2,}
  where the \eqn{G_j} are the random-effects covariance blocks.
}
\usage{
tolint_lmer(object, newdata = NULL, P = 0.90, conf = 0.95, side = 2,
            method = c("satterthwaite", "gpq", "beta-expectation"),
            B = 200, seed = NULL)
}
\arguments{
  \item{object}{a linear mixed model of class \code{"merMod"} fitted by
    \code{lme4::lmer()} (REML recommended).}
  \item{newdata}{a one-row data frame with the covariate values at which the
    interval is wanted. May be omitted only when the fixed part of the model
    is intercept-only and all random-effects terms are random intercepts.}
  \item{P}{content (for beta-content methods) or expected coverage (for
    \code{"beta-expectation"}), in (0, 1).}
  \item{conf}{confidence level for the beta-content methods, in (0, 1);
    ignored by \code{"beta-expectation"}.}
  \item{side}{1 for one-sided limits, 2 for a two-sided interval.}
  \item{method}{\code{"satterthwaite"}: beta-content interval that pairs a
    Satterthwaite-type effective degrees of freedom for
    \eqn{\hat\sigma_T^2} (its sampling variance estimated by parametric
    bootstrap) with the classical noncentral-t (one-sided) or Howe
    (two-sided) tolerance factors, in the spirit of Francq, Lin and Hoyer
    (2019); \code{"gpq"}: a parametric bootstrap analogue of the
    generalized pivotal quantity method; \code{"beta-expectation"}: a
    prediction-type interval whose coverage is \code{P} on average.}
  \item{B}{number of parametric bootstrap refits.}
  \item{seed}{optional integer seed for the bootstrap.}
}
\value{
  An object of class \code{"tolint"}; see \code{\link{tolint_oneway}}.
}
\details{
  All three methods refit the model to \code{B} parametric bootstrap
  samples generated from the fitted model, which can take a few seconds.
  Random-effects terms with covariates (random slopes) are supported as
  long as \code{newdata} provides the covariates appearing in those terms.
}
\references{
  Francq, B. G., Lin, D. and Hoyer, W. (2019). Confidence, prediction, and
  tolerance in linear mixed models. \emph{Statistics in Medicine}, 38,
  5603-5622.

  Krishnamoorthy, K. and Mathew, T. (2009). \emph{Statistical Tolerance
  Regions: Theory, Applications, and Computation}. Wiley.

  Hoffman, D. (2010). One-sided tolerance limits for balanced and
  unbalanced random effects models. \emph{Technometrics}, 52, 303-312.
}
\examples{
if (requireNamespace("lme4", quietly = TRUE)) {
  set.seed(2)
  k <- 10; n <- 4
  d <- data.frame(g = gl(k, n))
  d$y <- 5 + rep(rnorm(k, 0, 1.5), each = n) + rnorm(k * n)
  fit <- lme4::lmer(y ~ 1 + (1 | g), data = d)
  tolint_lmer(fit, side = 2, method = "satterthwaite", B = 100, seed = 1)
}
}
\keyword{htest}
