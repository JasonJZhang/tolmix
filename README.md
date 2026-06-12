# tolmix

Tolerance intervals for random effects and linear mixed-effects models in R.

A `(P, 1 − α)` **β-content tolerance interval** contains at least a
proportion `P` of a population with confidence `1 − α`. A
**β-expectation tolerance interval** contains a proportion `P` *on average*
(it is a prediction-type interval). This package computes both for the
marginal distribution of a future observation under variance-component
models — e.g. a future measurement from a *new* batch, lot, lab, or subject.

## Installation

```r
# install.packages("remotes")
remotes::install_github("JasonJZhang/tolmix")
```

`lme4` is needed only for `tolint_lmer()`.

## One-way random effects model

Model: `y_ij = mu + a_i + e_ij`, with `a_i ~ N(0, sigma_a^2)`,
`e_ij ~ N(0, sigma_e^2)`. The target is a future observation from a new
group: `Y ~ N(mu, sigma_a^2 + sigma_e^2)`. Balanced and unbalanced designs
are supported (unbalanced designs use the harmonic-mean adaptation based on
the unweighted mean of group means).

```r
library(tolmix)

set.seed(1)
g <- gl(12, 5)                               # 12 batches, 5 reps each
y <- 10 + rep(rnorm(12, 0, 2), each = 5) + rnorm(60)

# One-sided limits — Mee & Owen (1983), noncentral-t with Satterthwaite df
tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 1, method = "mee-owen")

# Two-sided interval — Hoffman & Kringle (2005) modified large sample (MLS)
tolint_oneway(y, g, P = 0.90, conf = 0.95, side = 2, method = "hoffman-kringle")

# Generalized pivotal quantities — Krishnamoorthy & Mathew (2004),
# one- or two-sided, Monte Carlo based
tolint_oneway(y, g, side = 2, method = "gpq", B = 50000)
```

## General linear mixed models (lme4)

For a model fitted with `lme4::lmer()`, the future observation is taken at
the fixed-effect covariates in `newdata` and from **new levels of every
grouping factor**, so its variance is the sum of all random-effect variance
contributions plus the residual variance. Random intercepts, random slopes,
multiple/crossed grouping factors, and fixed covariates are supported.

```r
library(lme4)
fit <- lmer(y ~ dose + (1 | batch) + (1 | analyst), data = mydata)

# Beta-content, Satterthwaite-type df (Francq, Lin & Hoyer 2019 spirit);
# Var(sigma_T^2_hat) is estimated by parametric bootstrap
tolint_lmer(fit, newdata = data.frame(dose = 50), side = 2,
            method = "satterthwaite", B = 500, seed = 1)

# Parametric-bootstrap pivotal (GPQ-analogue) method
tolint_lmer(fit, newdata = data.frame(dose = 50), side = 1,
            method = "gpq", B = 500, seed = 1)

# Beta-expectation (prediction-type) interval
tolint_lmer(fit, newdata = data.frame(dose = 50), side = 2,
            method = "beta-expectation", B = 500, seed = 1)
```

## Methods implemented

| Function | Method | Sides | Reference |
|---|---|---|---|
| `tolint_oneway` | `"mee-owen"` | 1 | Mee & Owen (1983), JASA 78, 901–905 |
| `tolint_oneway` | `"hoffman-kringle"` / `"mls"` | 2 | Hoffman & Kringle (2005), J. Biopharm. Stat. 15, 283–293 |
| `tolint_oneway` | `"gpq"` | 1 & 2 | Krishnamoorthy & Mathew (2004), Technometrics 46, 44–52 |
| `tolint_lmer` | `"satterthwaite"` | 1 & 2 | Francq, Lin & Hoyer (2019), Stat. Med. 38, 5603–5622 (spirit of) |
| `tolint_lmer` | `"gpq"` | 1 & 2 | parametric bootstrap analogue of Krishnamoorthy & Mathew (2009) |
| `tolint_lmer` | `"beta-expectation"` | 1 & 2 | Francq, Lin & Hoyer (2019) |

General reference: Krishnamoorthy, K. and Mathew, T. (2009),
*Statistical Tolerance Regions: Theory, Applications, and Computation*, Wiley.

## Notes

- For `side = 1`, both the lower and the upper one-sided limit are
  reported; each holds *marginally* at the stated confidence.
- The two-sided GPQ factor is found by root-finding on the Monte Carlo
  content-coverage criterion, following Krishnamoorthy & Mathew (2009).
- `tolint_lmer()` refits the model `B` times (parametric bootstrap), which
  takes a few seconds; increase `B` for more stable degrees of freedom.
