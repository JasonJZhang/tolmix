# Internal tolerance-factor utilities shared by the one-way and LMM methods.

# One-sided tolerance factor for X ~ N(mu, sigma^2) given an estimator
# mean ~ N(mu, sigma^2/ne) and an independent chi-square based variance
# estimate with f degrees of freedom (noncentral-t construction; see
# Krishnamoorthy and Mathew, 2009, Ch. 2 and Mee and Owen, 1983).
.k_onesided <- function(P, conf, ne, f) {
  delta <- qnorm(P) * sqrt(ne)
  suppressWarnings(qt(conf, df = f, ncp = delta)) / sqrt(ne)
}

# Two-sided beta-content factor via the Wald-Wolfowitz/Howe (1969)
# approximation: K = sqrt(f * chisq_P(1, ncp = 1/ne) / chisq_{1-conf}(f)).
.k_twosided <- function(P, conf, ne, f) {
  sqrt(f * qchisq(P, df = 1, ncp = 1 / ne) / qchisq(1 - conf, df = f))
}

# Modified large sample (Graybill-Wang) 1-conf upper confidence bound for a
# nonnegative linear combination sum(a * theta) of variance components, where
# est[i] estimates theta[i] with df[i] degrees of freedom.
.mls_upper <- function(a, est, df, conf) {
  H <- df / qchisq(1 - conf, df) - 1
  sum(a * est) + sqrt(sum((H * a * est)^2))
}

# Two-sided beta-content factor from a Monte Carlo set of generalized pivots
# (Gmu[i], Gsd[i]) for the mean and standard deviation of the sampled
# population: smallest K such that the interval center +/- K * scale has
# content >= P with (estimated) probability >= conf.
.k_gpq_twosided <- function(center, scale, Gmu, Gsd, P, conf) {
  cover <- function(K) {
    mean(pnorm((center + K * scale - Gmu) / Gsd) -
         pnorm((center - K * scale - Gmu) / Gsd) >= P) - conf
  }
  hi <- 2
  while (cover(hi) < 0 && hi < 1e4) hi <- hi * 2
  if (cover(hi) < 0) stop("failed to bracket the two-sided GPQ factor")
  uniroot(cover, lower = 0, upper = hi, tol = 1e-4)$root
}

# Constructor for the common return object.
.tolint <- function(model, method, side, P, conf, center, scale, K, limits,
                    estimates, info) {
  structure(list(model = model, method = method, side = side, P = P,
                 conf = conf, center = center, scale = scale, K = K,
                 limits = limits, estimates = estimates, info = info),
            class = "tolint")
}
