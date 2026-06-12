library(tolmix)

if (requireNamespace("lme4", quietly = TRUE)) {

  set.seed(2)
  k <- 12; n <- 5
  d <- data.frame(g = gl(k, n))
  d$y <- 5 + rep(rnorm(k, 0, 1.5), each = n) + rnorm(k * n, 0, 1)

  fit <- lme4::lmer(y ~ 1 + (1 | g), data = d)

  ## --- satterthwaite two-sided vs one-way Hoffman-Kringle ----------------
  tl <- tolint_lmer(fit, side = 2, method = "satterthwaite", B = 150, seed = 1)
  to <- tolint_oneway(d$y, d$g, side = 2, method = "hoffman-kringle")
  half_l <- unname(diff(tl$limits)) / 2
  half_o <- unname(diff(to$limits)) / 2
  # same target, similar widths (different df machinery, so allow 25%)
  stopifnot(abs(half_l - half_o) < 0.25 * half_o)
  stopifnot(abs(tl$center - to$center) < 0.2)
  # REML total variance close to ANOVA-type estimate on balanced data
  stopifnot(abs(tl$estimates["sigma.total2"] -
                to$estimates["sigma.total2"]) < 0.15 * to$estimates["sigma.total2"])

  ## --- one-sided satterthwaite vs Mee-Owen -------------------------------
  tl1 <- tolint_lmer(fit, side = 1, method = "satterthwaite", B = 150, seed = 1)
  to1 <- tolint_oneway(d$y, d$g, side = 1, method = "mee-owen")
  stopifnot(abs(tl1$limits["lower"] - to1$limits["lower"]) <
            0.3 * sqrt(to1$estimates["sigma.total2"]))

  ## --- gpq method runs and is ordered -------------------------------------
  tg <- tolint_lmer(fit, side = 2, method = "gpq", B = 150, seed = 3)
  stopifnot(tg$limits["lower"] < tg$center, tg$limits["upper"] > tg$center)

  ## --- beta-expectation is narrower than beta-content --------------------
  tb <- tolint_lmer(fit, side = 2, method = "beta-expectation", B = 150,
                    seed = 4)
  stopifnot(unname(diff(tb$limits)) < unname(diff(tl$limits)))

  ## --- model with a fixed covariate requires/uses newdata -----------------
  d$x <- rep(seq_len(n), k)
  d$y2 <- d$y + 0.5 * d$x
  fit2 <- lme4::lmer(y2 ~ x + (1 | g), data = d)
  res <- tryCatch(tolint_lmer(fit2, B = 60), error = function(e) "err")
  stopifnot(identical(res, "err"))
  t2 <- tolint_lmer(fit2, newdata = data.frame(x = 3), side = 2,
                    method = "satterthwaite", B = 100, seed = 5)
  stopifnot(abs(t2$center - (lme4::fixef(fit2)[1] + 3 * lme4::fixef(fit2)[2]))
            < 1e-8)

  cat("test-lmm.R: all checks passed\n")
} else {
  cat("test-lmm.R: lme4 not available, skipped\n")
}
