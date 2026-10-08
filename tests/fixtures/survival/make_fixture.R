# The R oracle for ctx.clinical.survival and its page twin (GS plan, test 3): made-up cohorts
# drawn here, and what R's survival package computes on them, written to fixture.json.
# R is not in pixi; this script and its output are committed, and the output records the R
# and survival versions. Rerun: Rscript tests/fixtures/survival/make_fixture.R
suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
})
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
set.seed(20261008)

# Cohort A: 160 patients, times in days rounded to 30-day steps so events tie, two
# covariates (a z-like expression and a stage 1 to 4), a two- and a three-level group.
n <- 160
x1 <- round(rnorm(n), 4)
x2 <- sample(1:4, n, replace = TRUE, prob = c(0.3, 0.4, 0.2, 0.1))
hazard <- 0.0008 * exp(0.4 * x1 + 0.35 * (x2 - 1))
event_time <- rexp(n, hazard)
censor_time <- runif(n, 200, 3000)
time <- pmax(30, round(pmin(event_time, censor_time) / 30) * 30)
event <- as.integer(event_time <= censor_time)
g2 <- ifelse(x1 >= median(x1), "high", "low")
g3 <- cut(x1, quantile(x1, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE, labels = c("t1", "t2", "t3"))
a <- data.frame(time, event, x1, x2, g2, g3 = as.character(g3))
# Three time-zero rows and a time tied between an event and a censored patient.
a$time[1:3] <- 0
a$event[1:3] <- c(1, 0, 0)

km <- survfit(Surv(time, event) ~ 1, data = a, conf.type = "log-log")
km_table <- summary(km, censored = TRUE)
med <- quantile(km, probs = 0.5)
horizon <- 1000
rm <- summary(km, rmean = horizon)$table
reverse <- survfit(Surv(time, 1 - event) ~ 1, data = a)
km_by <- survfit(Surv(time, event) ~ g2, data = a, conf.type = "log-log")
med_by <- quantile(km_by, probs = 0.5)

fit_cox <- function(formula, ties) {
  fit <- coxph(formula, data = a, ties = ties)
  zph <- cox.zph(fit, transform = "rank")
  list(coef = unname(coef(fit)), se = unname(sqrt(diag(vcov(fit)))),
       loglik = unname(fit$loglik), concordance = unname(fit$concordance["concordance"]),
       zph_chisq = unname(zph$table[, "chisq"]), zph_p = unname(zph$table[, "p"]),
       zph_terms = rownames(zph$table))
}
lr2 <- survdiff(Surv(time, event) ~ g2, data = a)
lr3 <- survdiff(Surv(time, event) ~ g3, data = a)

# Cohort B: one patient with a non-zero value, censored, as ENSG00000248769 was on the GS
# plan's lists (gene_selection_plan.md): the likelihood rises without bound as beta falls.
b <- data.frame(time = c(seq(30, 1800, by = 30)), event = rep(c(1, 0, 0), 20), x = 0)
b$x[2] <- 1.5
monotone_warning <- ""
fit_b <- withCallingHandlers(
  coxph(Surv(time, event) ~ x, data = b),
  warning = function(w) {
    monotone_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)

out <- list(
  r_version = R.version.string, survival_version = as.character(packageVersion("survival")),
  seed = 20261008, horizon = horizon, level = 0.95,
  a = a,
  km = list(time = km_table$time, n_risk = km_table$n.risk, n_event = km_table$n.event,
            n_censor = km_table$n.censor, surv = km_table$surv, lower = km_table$lower,
            upper = km_table$upper, median = unname(med$quantile), median_lower = unname(med$lower),
            median_upper = unname(med$upper), rmean = unname(rm["rmean"]),
            reverse_median = unname(quantile(reverse, probs = 0.5)$quantile)),
  km_by_g2 = list(strata = names(km_by$strata), median = unname(med_by$quantile),
                  median_lower = unname(med_by$lower), median_upper = unname(med_by$upper)),
  logrank_g2 = list(chisq = unname(lr2$chisq), df = length(lr2$n) - 1,
                    p = pchisq(lr2$chisq, length(lr2$n) - 1, lower.tail = FALSE)),
  logrank_g3 = list(chisq = unname(lr3$chisq), df = length(lr3$n) - 1,
                    p = pchisq(lr3$chisq, length(lr3$n) - 1, lower.tail = FALSE)),
  cox_x1_efron = fit_cox(Surv(time, event) ~ x1, "efron"),
  cox_x1x2_efron = fit_cox(Surv(time, event) ~ x1 + x2, "efron"),
  cox_x1x2_breslow = fit_cox(Surv(time, event) ~ x1 + x2, "breslow"),
  cox_x1_strata = {
    fit <- coxph(Surv(time, event) ~ x1 + strata(g2), data = a)
    zph <- cox.zph(fit, transform = "rank")
    list(coef = unname(coef(fit)), se = unname(sqrt(diag(vcov(fit)))),
         zph_chisq = unname(zph$table[, "chisq"]), zph_p = unname(zph$table[, "p"]))
  },
  b = b,
  monotone = list(coef = unname(coef(fit_b)), se = unname(sqrt(diag(vcov(fit_b)))),
                  warning = monotone_warning)
)
writeLines(toJSON(out, digits = NA, auto_unbox = TRUE, pretty = TRUE),
           file.path(here, "fixture.json"))
