# The R oracle for ctx.clinical.survival.cumulative_incidence (M10.11): a made-up cohort with
# two competing causes and ties, and the Aalen-Johansen estimate R's survival package gives,
# survfit on a factor status, written to cif.json with the R and survival versions.
# R is not in pixi; this script and its output are committed.
# Rerun: Rscript tests/fixtures/survival/make_cif.R
suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
})
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
set.seed(20261009)
n <- 120
time <- pmax(30, round(rexp(n, 1 / 900) / 30) * 30)
cause <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
fit <- survfit(Surv(time, factor(cause, 0:2, c("censor", "one", "two"))) ~ 1)
out <- list(
  r = R.version.string, survival = as.character(packageVersion("survival")),
  time = time, cause = cause, times = fit$time,
  one = unname(fit$pstate[, which(fit$states == "one")]),
  two = unname(fit$pstate[, which(fit$states == "two")])
)
writeLines(toJSON(out, digits = NA, auto_unbox = TRUE, pretty = TRUE), file.path(here, "cif.json"))
