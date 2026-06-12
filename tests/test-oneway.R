library(tolmix)

set.seed(1)
k <- 15; n <- 6
g <- gl(k, n)
y <- 10 + rep(rnorm(k, 0, 2), each = n) + rnorm(k * n, 0, 1)

zp <- qnorm(0.90)

## --- Mee-Owen one-sided -----------------------------------------------
t1 <- tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 1, method = "mee-owen")
stopifnot(t1$limits["lower"] < t1$center, t1$limits["upper"] > t1$center)
# the tolerance factor must exceed the naive normal quantile
stopifnot(t1$K > zp)
# variance estimates are coherent
stopifnot(abs(t1$estimates["sigma.total2"] -
              (t1$estimates["sigma.between2"] + t1$estimates["sigma.within2"])) < 1e-8)

## --- Hoffman-Kringle two-sided ----------------------------------------
t2 <- tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 2,
                    method = "hoffman-kringle")
stopifnot(t2$K > qnorm(0.95))           # wider than the naive z interval
stopifnot(t2$limits["lower"] < t1$limits["lower"] + 5)  # sanity, finite
stopifnot(is.finite(t2$limits["lower"]), is.finite(t2$limits["upper"]))
# "mls" alias gives the same result
t2b <- tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 2, method = "mls")
stopifnot(identical(t2$limits, t2b$limits))

## --- GPQ vs Mee-Owen (one-sided) should roughly agree ------------------
set.seed(42)
t3 <- tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 1, method = "gpq",
                    B = 50000)
sdy <- sqrt(t1$estimates["sigma.total2"])
stopifnot(abs(t3$limits["lower"] - t1$limits["lower"]) < 0.25 * sdy)
stopifnot(abs(t3$limits["upper"] - t1$limits["upper"]) < 0.25 * sdy)

## --- GPQ two-sided vs Hoffman-Kringle should roughly agree -------------
set.seed(43)
t4 <- tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 2, method = "gpq",
                    B = 50000)
stopifnot(abs(t4$K - t2$K) < 0.35)

## --- large samples: K approaches the normal quantile -------------------
set.seed(7)
k2 <- 400; n2 <- 10
g2 <- gl(k2, n2)
y2 <- rep(rnorm(k2, 0, 1), each = n2) + rnorm(k2 * n2, 0, 1)
t5 <- tolint_oneway(y2, g2, P = 0.90, conf = 0.95, side = 1,
                    method = "mee-owen")
stopifnot(t5$K < zp * 1.12, t5$K > zp)

## --- unbalanced design runs and is sane ---------------------------------
set.seed(8)
ni <- sample(2:8, 12, replace = TRUE)
gu <- factor(rep(seq_along(ni), ni))
yu <- 3 + rep(rnorm(12, 0, 1.5), ni) + rnorm(sum(ni))
u1 <- tolint_oneway(yu, gu, side = 1, method = "mee-owen")
u2 <- tolint_oneway(yu, gu, side = 2, method = "hoffman-kringle")
stopifnot(u1$K > zp, u2$K > qnorm(0.95),
          u2$limits["upper"] > u2$limits["lower"])

## --- argument validation -------------------------------------------------
res <- tryCatch(tolint_oneway(y, g, side = 2, method = "mee-owen"),
                error = function(e) "err")
stopifnot(identical(res, "err"))
res <- tryCatch(tolint_oneway(y, g, side = 1, method = "hoffman-kringle"),
                error = function(e) "err")
stopifnot(identical(res, "err"))

cat("test-oneway.R: all checks passed\n")
