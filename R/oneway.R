# Tolerance intervals for the one-way random effects model
#   y_ij = mu + a_i + e_ij,  a_i ~ N(0, sigma_a^2),  e_ij ~ N(0, sigma_e^2),
# for the marginal distribution of a future observation from a new group,
#   Y ~ N(mu, sigma_a^2 + sigma_e^2).
#
# Balanced and unbalanced designs are supported; unbalanced designs use the
# harmonic-mean adaptation of Krishnamoorthy and Mathew (2004), based on the
# unweighted mean of the group means.

tolint_oneway <- function(y, group, P = 0.90, conf = 0.95, side = 2,
                          method = c("auto", "mee-owen", "hoffman-kringle",
                                     "gpq", "mls"),
                          B = 10000) {
  method <- match.arg(method)
  if (method == "mls") method <- "hoffman-kringle"
  if (!is.numeric(y)) stop("'y' must be numeric")
  if (length(y) != length(group)) stop("'y' and 'group' lengths differ")
  if (!(side %in% c(1, 2))) stop("'side' must be 1 or 2")
  if (P <= 0 || P >= 1 || conf <= 0 || conf >= 1)
    stop("'P' and 'conf' must be in (0, 1)")

  ok <- is.finite(y) & !is.na(group)
  y <- y[ok]
  group <- factor(group[ok])
  k <- nlevels(group)
  N <- length(y)
  if (k < 2) stop("need at least 2 groups")
  if (N - k < 1) stop("need at least one group with replication")

  if (method == "auto")
    method <- if (side == 1) "mee-owen" else "hoffman-kringle"
  if (method == "mee-owen" && side != 1)
    stop("method 'mee-owen' is one-sided only; use 'hoffman-kringle' or 'gpq'")
  if (method == "hoffman-kringle" && side != 2)
    stop("method 'hoffman-kringle' is two-sided only; use 'mee-owen' or 'gpq'")

  ni <- as.vector(table(group))
  ntilde <- k / sum(1 / ni)              # harmonic mean group size
  ybari <- tapply(y, group, mean)
  ybar <- mean(ybari)                    # unweighted mean of group means
  S1sq <- sum((ybari - ybar)^2) / (k - 1)  # estimates sigma_a^2 + sigma_e^2/ntilde
  f1 <- k - 1
  sse <- sum((y - ybari[group])^2)
  f2 <- N - k
  S2sq <- sse / f2                       # estimates sigma_e^2
  sigy2 <- S1sq + (1 - 1 / ntilde) * S2sq  # estimates sigma_a^2 + sigma_e^2
  siga2 <- max(S1sq - S2sq / ntilde, 0)

  zp <- qnorm(P)
  scale <- sqrt(sigy2)

  if (method == "mee-owen") {
    # Effective sample size for ybar relative to sigma_y^2, Satterthwaite df.
    ne <- k * sigy2 / S1sq
    f <- sigy2^2 / (S1sq^2 / f1 + ((1 - 1 / ntilde) * S2sq)^2 / f2)
    K <- .k_onesided(P, conf, ne, f)
    limits <- c(lower = ybar - K * scale, upper = ybar + K * scale)
  } else if (method == "hoffman-kringle") {
    # MLS upper bound on gamma = Var(ybar) + sigma_y^2
    #   = (1 + 1/k) * (sigma_a^2 + sigma_e^2/ntilde) + (1 - 1/ntilde) * sigma_e^2.
    a <- c(1 + 1 / k, 1 - 1 / ntilde)
    gamU <- .mls_upper(a, c(S1sq, S2sq), c(f1, f2), conf)
    half <- qnorm((1 + P) / 2) * sqrt(gamU)
    K <- half / scale
    limits <- c(lower = ybar - half, upper = ybar + half)
  } else { # gpq
    Z <- rnorm(B)
    U1 <- rchisq(B, f1)
    U2 <- rchisq(B, f2)
    Gs1 <- f1 * S1sq / U1                # pivot for sigma_a^2 + sigma_e^2/ntilde
    Gse <- f2 * S2sq / U2                # pivot for sigma_e^2
    Gsy <- pmax(Gs1 - Gse / ntilde, 0) + Gse
    Gmu <- ybar - Z * sqrt(Gs1 / k)
    if (side == 1) {
      K <- NA_real_
      limits <- c(lower = unname(quantile(Gmu - zp * sqrt(Gsy), 1 - conf)),
                  upper = unname(quantile(Gmu + zp * sqrt(Gsy), conf)))
    } else {
      K <- .k_gpq_twosided(ybar, scale, Gmu, sqrt(Gsy), P, conf)
      limits <- c(lower = ybar - K * scale, upper = ybar + K * scale)
    }
  }

  est <- c(mean = unname(ybar), sigma.between2 = siga2,
           sigma.within2 = S2sq, sigma.total2 = sigy2)
  info <- sprintf("k = %d groups, N = %d observations, %s design",
                  k, N, if (length(unique(ni)) == 1) "balanced" else "unbalanced")
  .tolint("one-way random effects model", method, side, P, conf,
          unname(ybar), scale, K, limits, est, info)
}
