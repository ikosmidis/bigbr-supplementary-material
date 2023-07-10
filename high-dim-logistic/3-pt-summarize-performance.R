library("dplyr")
options(dplyr.summarise.inform = FALSE)

if (interactive()) {
    b_setting <- "a"
    nobs <- 2000
    rhosq <- 0.0
    base_path <- "~/Repositories/bigbr-supplementary-material/high-dim-logistic/"
}
results_path <- file.path(base_path, "results")
figures_path <- file.path(base_path, "figures")
source(file.path(base_path, "functions.R"))

plot_type <- if (b_setting == "a") "estimate_vs_truth" else "estimate_and_truth"


base_name <- paste0("mJPL-", rhosq, "-n-", nobs, "-setting-", b_setting)
## Get phase transfition curve for rhosq
load(file.path(results_path, paste0("PT-n-", nobs, "-setting-", b_setting, "-rhosq-", rhosq, ".rda")))
## Get estimates for rhosq
files <- dir(results_path,
             pattern = paste0("estimates-n-", nobs, "-setting-", b_setting, "-rhosq-", rhosq, "-"),
             full.names = TRUE)
c_perf <- NULL
for (f in files) {
    load(f)
    c_perf <- rbind(c_perf, attr(results, "performance"))
}

perf <- c_perf |> subset(method == "mJPL") |>
    group_by(kappa, gamma, p) |>
    summarize(mean_elapsed = mean(elapsed),
              min_iter = min(iter),
              avg_iter = mean(iter),
              max_iter = max(iter)) |> data.frame()



print(perf)

## library("memisc")
## toLatex(perf, digits = 2)
