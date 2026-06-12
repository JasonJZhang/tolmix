\name{tolint_oneway}
\alias{tolint_oneway}
\title{Tolerance intervals for the one-way random effects model}
\description{
  Computes one-sided or two-sided \eqn{(P, 1-\alpha)} beta-content tolerance
  intervals for the marginal distribution of a future observation
  \eqn{Y \sim N(\mu, \sigma_a^2 + \sigma_e^2)} under the one-way random
  effects model \eqn{y_{ij} = \mu + a_i + e_{ij}}, with
  \eqn{a_i \sim N(0, \sigma_a^2)} and \eqn{e_{ij} \sim N(0, \sigma_e^2)}.
  Balanced and unbalanced designs are supported; unbalanced designs use the
  harmonic-mean adaptation based on the unweighted mean of the group means
  (Krishnamoorthy and Mathew, 2004).
}
\usage{
tolint_oneway(y, group, P = 0.90, conf = 0.95, side = 2,
              method = c("auto", "mee-owen", "hoffman-kringle", "gpq", "mls"),
              B = 10000)
}
\arguments{
  \item{y}{numeric response vector.}
  \item{group}{vector or factor of the same length as \code{y} identifying
    the random groups (batches, lots, subjects, ...).}
  \item{P}{content of the tolerance interval (proportion of the population
    to be covered), in (0, 1).}
  \item{conf}{confidence level \eqn{1-\alpha}, in (0, 1).}
  \item{side}{1 for one-sided limits, 2 for a two-sided interval.}
  \item{method}{\code{"mee-owen"}: one-sided limits of Mee and Owen (1983)
    based on a noncentral-t factor with Satterthwaite degrees of freedom;
    \code{"hoffman-kringle"} (alias \code{"mls"}): two-sided interval of
    Hoffman and Kringle (2005) based on a modified large sample upper bound;
    \code{"gpq"}: generalized pivotal quantity method of Krishnamoorthy and
    Mathew (2004), available for both \code{side = 1} and \code{side = 2}
    (Monte Carlo). \code{"auto"} picks \code{"mee-owen"} for \code{side = 1}
    and \code{"hoffman-kringle"} for \code{side = 2}.}
  \item{B}{number of Monte Carlo samples for \code{method = "gpq"}.}
}
\value{
  An object of class \code{"tolint"}: a list with the limits, the tolerance
  factor (when defined), variance component estimates and settings. For
  \code{side = 1} both the lower and the upper one-sided limit are reported;
  each holds marginally at confidence \code{conf}.
}
\references{
  Mee, R. W. and Owen, D. B. (1983). Improved factors for one-sided
  tolerance limits for balanced one-way ANOVA random model.
  \emph{Journal of the American Statistical Association}, 78, 901-905.

  Krishnamoorthy, K. and Mathew, T. (2004). One-sided tolerance limits in
  balanced and unbalanced one-way random models based on generalized
  confidence intervals. \emph{Technometrics}, 46, 44-52.

  Hoffman, D. and Kringle, R. (2005). Two-sided tolerance intervals for
  balanced and unbalanced random effects models.
  \emph{Journal of Biopharmaceutical Statistics}, 15, 283-293.

  Krishnamoorthy, K. and Mathew, T. (2009). \emph{Statistical Tolerance
  Regions: Theory, Applications, and Computation}. Wiley.
}
\examples{
set.seed(1)
k <- 12; n <- 5
g <- gl(k, n)
y <- 10 + rep(rnorm(k, 0, 2), each = n) + rnorm(k * n)

# one-sided limits (Mee-Owen)
tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 1)

# two-sided interval (Hoffman-Kringle MLS)
tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 2)

# generalized pivotal quantity method
tolint_oneway(y, g, side = 2, method = "gpq", B = 5000)
}
\keyword{htest}
