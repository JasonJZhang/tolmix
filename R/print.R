print.tolint <- function(x, digits = 4, ...) {
  labels <- c("mee-owen" = "Mee-Owen (1983) one-sided noncentral-t",
              "hoffman-kringle" = "Hoffman-Kringle (2005) modified large sample",
              "gpq" = "generalized pivotal quantity (Monte Carlo)",
              "satterthwaite" = "Satterthwaite-type df with Howe/noncentral-t factor",
              "beta-expectation" = "beta-expectation (prediction-type)")
  lab <- labels[x$method]
  if (is.na(lab)) lab <- x$method
  cat(sprintf("%d-sided %s tolerance interval\n", x$side,
              if (x$method == "beta-expectation") "beta-expectation"
              else "beta-content"))
  cat("Model:     ", x$model, "\n", sep = "")
  cat("Method:    ", lab, "\n", sep = "")
  cat(sprintf("Content P = %.3f, Confidence = %.3f\n", x$P, x$conf))
  cat(x$info, "\n", sep = "")
  cat("\nEstimates:\n")
  print(round(x$estimates, digits))
  if (!is.na(x$K))
    cat(sprintf("\nTolerance factor K = %.*f (relative to sigma.total)\n",
                digits, x$K))
  cat("\n")
  if (x$side == 2) {
    cat(sprintf("Two-sided tolerance interval: [%.*f, %.*f]\n",
                digits, x$limits["lower"], digits, x$limits["upper"]))
  } else {
    cat(sprintf("One-sided lower tolerance limit: %.*f\n",
                digits, x$limits["lower"]))
    cat(sprintf("One-sided upper tolerance limit: %.*f\n",
                digits, x$limits["upper"]))
    cat("(each one-sided limit holds marginally at the stated confidence)\n")
  }
  invisible(x)
}
