## Rscript --no-init-file -e 'nobs <- 500; beta_star_setting <- "a"; ncores <- 5; source("1-pt-simulate.R")'
## Rscript --no-init-file -e 'nobs <- 500; beta_star_setting <- "b"; ncores <- 5; source("1-pt-simulate.R")'
## Rscript --no-init-file -e 'nobs <- 2000; beta_star_setting <- "a"; ncores <- 5; source("1-pt-simulate.R")'
## Rscript --no-init-file -e 'nobs <- 2000; beta_star_setting <- "b"; ncores <- 5; source("1-pt-simulate.R")'

devtools::load_all("~/Repositories/biglm")
library("brglm2")
library("parallel")

code_path <- "~/Repositories/bigbr-supplementary-material/high-dim-logistic/"
results_path <- file.path(code_path, "results")
source(file.path(code_path, "functions.R"))

## beta_star_setting and nobs are set by the following lines only if R
## is used interactively
if (interactive()) {
    beta_star_setting <- "a"
    nobs <- 2000
    ncores <- 5
}

ns <- 200000
repetitions <- 5
maxit <- 250
tolerance <- 1e-03
chunksize <- 1000
rhosq_grid <- c(0.0, 0.75, 0.25, 0.5, 0.8, 0.9)

base_settings <- data.frame(kappa = c(0.01, 0.01, 0.01,
                                      0.05, 0.05, 0.05,
                                      0.15, 0.15, 0.15,
                                      0.22, 0.22,
                                      0.25, 0.25, 0.25,
                                      0.30, 0.30,
                                      0.35,
                                      0.35, 0.35, 0.35,
                                      0.40, 0.40,
                                      0.45, 0.45, 0.45,
                                      0.50, 0.50,
                                      0.55, 0.55, 0.55
                                      ),
                            gamma = c(1, 8, 15,
                                      4.5, 11.5, 18.5,
                                      1, 8, 15,
                                      8, 15,
                                      4.5, 11.5, 18.5,
                                      8, 15,
                                      1,
                                      4.5, 11.5, 18.5,
                                      8, 15,
                                      4.5, 11.5, 18.5,
                                      8, 15,
                                      4.5, 11.5, 18.5
                                      ),
                            boundary = FALSE)


all_settings <- base_settings
all_settings$n <- nobs
all_settings$p <- ceiling(all_settings$n * all_settings$kappa)
all_settings$repetitions <- repetitions
n_settings <- nrow(all_settings)
all_settings$maxit <- maxit
all_settings$tolerance <- tolerance
all_settings$chunksize <- chunksize

set.seed(111)

for (rhosq in rhosq_grid) {

    ## Compute phase transition for the particular value fo rhosq
    cat("ρ^2 = ", rhosq, "\n")
    xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
    ga <- c(0.001, 0.01, seq(0, 20, length = 200))
    pt <- mclapply(1:length(ga), function(k) {
        d <- compute_PT(beta0 = sqrt(rhosq) * ga[k], ga[k], ncores = 1, XZU =  xzu)
        d$rhosq <- rhosq
        d
    }, mc.cores = ncores)
    pt <- do.call("rbind", pt)

    ## Save phase transition curve
    out_path <- file.path(results_path, paste0("PT-n-", nobs, "-setting-", beta_star_setting, "-rhosq-", rhosq, ".rda"))
    save(pt, file = out_path)

    ## Determine whether ML estimates exist or not and set up and add
    ## to settings
    all_settings$rhosq <- rhosq
    all_settings$mle_exists <- all_settings$kappa <= approx(pt[2:1], xout = all_settings$gamma)$y
    all_settings$mle_exists[all_settings$boundary] <- FALSE
    ## Make sure that unique seeds are used
    dup_check <- TRUE
    while (dup_check) {
        all_settings$seed <- round(runif(n_settings) * 1000000)
        dup_check <- length(unique(all_settings$seed)) != n_settings
    }

    ## Compute estimates
    for (wh in seq.int(n_settings)) {
        results <- get_results(all_settings[wh, ],
                               beta0 = sqrt(rhosq) * all_settings[wh, "gamma"],
                               nz_perc = 0.2,
                               beta_star_setting = beta_star_setting,
                               ncores = ncores,
                               verbose = 0)
        out_path <- file.path(results_path, paste0("estimates-n-", nobs, "-setting-",
                                                   beta_star_setting, "-rhosq-", rhosq,
                                                   "-kappa-", round(all_settings[wh, "kappa"], 4),
                                                   "-gamma-", all_settings[wh, "gamma"],
                                                   ".rda"))
        save(results, file = out_path)
    }

}
