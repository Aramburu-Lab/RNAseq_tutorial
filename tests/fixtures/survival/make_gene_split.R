# The R oracle for ctx.clinical.gene.split (M10.G): fixture.json's cohort A split at the median
# of four made-up genes, at or above the median against below, and what R's survival package
# computes on each half: the Kaplan-Meier steps (conf.type "log-log"), the number at risk at
# fixed times, and survdiff's log-rank. Written to gene_split.json with the R and survival
# versions. R is not in pixi; rerun: Rscript tests/fixtures/survival/make_gene_split.R
suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
})
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
a <- fromJSON(file.path(here, "fixture.json"))$a

# x1 is continuous; x2, a stage 1 to 4, ties at its median; tied rounds x1, so many patients
# sit on the median; floor holds more than half its patients at 0, its median, so no patient
# is below it and there is no test.
genes <- list(x1 = a$x1, x2 = a$x2, tied = round(a$x1), floor = pmax(round(a$x1), 0))
ticks <- seq(0, max(a$time), by = 360)

half <- function(rows) {
  if (!nrow(rows)) return(NULL)
  fit <- survfit(Surv(time, event) ~ 1, data = rows, conf.type = "log-log")
  s <- summary(fit, censored = TRUE)
  list(time = s$time, n_risk = s$n.risk, n_event = s$n.event, n_censor = s$n.censor,
       surv = s$surv, lower = s$lower, upper = s$upper,
       risk = summary(fit, times = ticks, extend = TRUE)$n.risk)
}

out <- list(r_version = R.version.string,
            survival_version = as.character(packageVersion("survival")), ticks = ticks)
for (name in names(genes)) {
  x <- genes[[name]]
  upper <- as.integer(x >= median(x))
  rows <- data.frame(time = a$time, event = a$event, upper = upper)
  test <- NULL
  if (length(unique(upper)) == 2) {
    lr <- survdiff(Surv(time, event) ~ upper, data = rows)
    test <- list(chisq = unname(lr$chisq), p = pchisq(lr$chisq, 1, lower.tail = FALSE))
  }
  out[[name]] <- list(values = x, median = median(x), upper = upper,
                      below = half(rows[rows$upper == 0, ]), above = half(rows[rows$upper == 1, ]),
                      logrank = test)
}
write_json(out, file.path(here, "gene_split.json"), digits = NA, auto_unbox = TRUE, null = "null",
           na = "string")
