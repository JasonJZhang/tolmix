# Tolerance intervals for general linear mixed models fitted with
# lme4::lmer(). The target distribution is that of a future observation at
# fixed-effect covariates 'newdata', drawn from NEW levels of every grouping
# factor:
#   Y ~ N(x0' beta, sigma_T^2),  sigma_T^2 = sum_j z0j' G_j z0j + sigma_e^2.
#
# Methods:
#   "satterthwaite"   - beta-content interval using a Satterthwaite-type
#                       effective df for sigma_T^2 (variance of the estimate
#                       obtained by parametric bootstrap), combined with the
#                       Mee-Owen one-sided / Howe two-sided factors; in the
#                       spirit of Francq, Lin and Hoyer (2019).
#   "gpq"             - parametric bootstrap analogue of the generalized
#                       pivotal quantity method of Krishnamoorthy and Mathew.
#   "beta-expectation"- beta-expectation interval (mean coverage P), i.e. a
#                       prediction-type interval x0'beta_hat +/-
#                       t * sqrt(Var(x0'beta_hat) + sigma_T^2).

tolint_lmer <- function(object, newdata = NULL, P = 0.90, conf = 0.95,
                        side = 2,
                        method = c("satterthwaite", "gpq", "beta-expectation"),
                        B = 200, seed = NULL) {
  method <- match.arg(method)
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("package 'lme4' is required for tolint_lmer()")
  if (!inherits(object, "merMod") || !lme4::isLMM(object))
    stop("'object' must be a linear mixed model fitted with lme4::lmer()")
  if (!(side %in% c(1, 2))) stop("'side' must be 1 or 2")
  if (P <= 0 || P >= 1 || conf <= 0 || conf >= 1)
    stop("'P' and 'conf' must be in (0, 1)")
  if (B < 50) stop("'B' should be at least 50")

  pars <- .lmm_pars(object, newdata)
  mu <- pars$mu
  varmu <- pars$varmu
  sigT2 <- pars$sigT2
  scale <- sqrt(sigT2)
  zp <- qnorm(P)

  # Parametric bootstrap of (mu_hat, Var(mu_hat), sigma_T^2_hat).
  if (!is.null(seed)) set.seed(seed)
  sims <- simulate(object, nsim = B)
  boot <- vapply(sims, function(yb) {
    fit <- tryCatch(suppressWarnings(suppressMessages(lme4::refit(object, yb))),
                    error = function(e) NULL)
    if (is.null(fit)) return(c(NA_real_, NA_real_, NA_real_))
    p <- tryCatch(.lmm_pars(fit, newdata), error = function(e) NULL)
    if (is.null(p)) return(c(NA_real_, NA_real_, NA_real_))
    c(p$mu, p$varmu, p$sigT2)
  }, numeric(3))
  boot <- boot[, complete.cases(t(boot)), drop = FALSE]
  if (ncol(boot) < 0.5 * B)
    warning("more than half of the bootstrap refits failed")
  if (ncol(boot) < 30)
    stop("too few successful bootstrap refits to proceed")

  if (method == "satterthwaite") {
    vT <- var(boot[3, ])
    nu <- 2 * sigT2^2 / vT
    ne <- sigT2 / varmu
    K <- if (side == 1) .k_onesided(P, conf, ne, nu)
         else .k_twosided(P, conf, ne, nu)
    limits <- c(lower = mu - K * scale, upper = mu + K * scale)
    extra <- c(df.satt = nu, n.eff = ne)
  } else if (method == "gpq") {
    Gmu <- 2 * mu - boot[1, ]                 # basic bootstrap pivot for mu
    GsT2 <- sigT2^2 / boot[3, ]               # ratio pivot for sigma_T^2
    if (side == 1) {
      K <- NA_real_
      limits <- c(lower = unname(quantile(Gmu - zp * sqrt(GsT2), 1 - conf)),
                  upper = unname(quantile(Gmu + zp * sqrt(GsT2), conf)))
    } else {
      K <- .k_gpq_twosided(mu, scale, Gmu, sqrt(GsT2), P, conf)
      limits <- c(lower = mu - K * scale, upper = mu + K * scale)
    }
    extra <- c(B.used = ncol(boot))
  } else { # beta-expectation
    s2 <- varmu + sigT2
    vs2 <- var(boot[2, ] + boot[3, ])
    nu <- 2 * s2^2 / vs2
    tq <- if (side == 1) qt(P, nu) else qt((1 + P) / 2, nu)
    K <- tq * sqrt(s2) / scale
    limits <- c(lower = mu - tq * sqrt(s2), upper = mu + tq * sqrt(s2))
    extra <- c(df.satt = nu)
  }

  est <- c(mean = mu, var.mean = varmu, sigma.total2 = sigT2,
           sigma.resid2 = pars$sige2, extra)
  info <- sprintf("lmer fit, %d obs; bootstrap B = %d (%d successful)",
                  length(stats::fitted(object)), B, ncol(boot))
  .tolint("linear mixed-effects model (lme4)", method, side, P, conf,
          mu, scale, K, limits, est, info)
}

# Extract (mu, Var(mu), sigma_T^2) at 'newdata' from an lmer fit.
.lmm_pars <- function(object, newdata) {
  beta <- lme4::fixef(object)
  tt <- delete.response(terms(object))
  if (is.null(newdata)) {
    if (length(beta) != 1L || names(beta)[1L] != "(Intercept)")
      stop("'newdata' (a one-row data frame) is required unless the fixed ",
           "part of the model is intercept-only")
    X0 <- matrix(1, 1, 1)
  } else {
    if (NROW(newdata) != 1L)
      stop("'newdata' must have exactly one row")
    mf <- model.frame(tt, newdata,
                      xlev = .getXlevels(tt, model.frame(object)))
    X0 <- model.matrix(tt, mf,
                       contrasts.arg = attr(model.matrix(object), "contrasts"))
  }
  if (ncol(X0) != length(beta))
    stop("could not align 'newdata' with the fixed-effects design")
  mu <- drop(X0 %*% beta)
  varmu <- drop(X0 %*% as.matrix(vcov(object)) %*% t(X0))

  vc <- lme4::VarCorr(object)
  bars <- if (requireNamespace("reformulas", quietly = TRUE))
    reformulas::findbars(formula(object))
  else
    suppressWarnings(lme4::findbars(formula(object)))
  if (length(vc) != length(bars))
    stop("unsupported random-effects structure")
  vre <- 0
  for (i in seq_along(bars)) {
    G <- as.matrix(vc[[i]])
    lhs <- paste(deparse(bars[[i]][[2]]), collapse = "")
    if (lhs == "1") {
      z0 <- matrix(1, 1, 1)
    } else {
      if (is.null(newdata))
        stop("'newdata' is required for random-effects terms with covariates")
      z0 <- model.matrix(as.formula(paste("~", lhs)), newdata)
      if (!all(colnames(G) %in% colnames(z0)))
        stop("could not align 'newdata' with random-effects term '", lhs, "'")
      z0 <- z0[, colnames(G), drop = FALSE]
    }
    vre <- vre + drop(z0 %*% G %*% t(z0))
  }
  sige2 <- sigma(object)^2
  list(mu = mu, varmu = varmu, sige2 = sige2, sigT2 = vre + sige2)
}
