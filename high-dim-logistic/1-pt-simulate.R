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



if (FALSE) {
    ## Setting for Figure 2b on p 11 of the supplementary information appendix of
    ##
    ## Sur P, and Candès EJ (2019). A Modern Maximum-Likelihood Theory for
    ## High-Dimensional Logistic Regression.  Proceedings of the National
    ## Academy of Sciences 116 (29):
    ## 14516–25. https://doi.org/10.1073/pnas.1810420116.
    ##
    ## The supplementary information appendix has been downloaded from
    ## https://www.pnas.org/content/pnas/suppl/2019/06/29/1810420116.DCSupplemental/pnas.1810420116.sapp.pdf

    set.seed(111)
    n <- 2000
    p <- 400
    ## beta <- c(rep(10, p/8), rep(-10, p/8), rep(0, 3*p/4))
    beta <- c(rep(10, p/8), rep(5, p/8), rep(0, p/4), rep(-2, p/4), rep(-15, p/4))
    x <- matrix(rnorm(n * p, 0, sqrt(1/n)), n, p)
    probs <- plogis(x %*% beta)
    y <- rbinom(n, 1, probs)
    form <- formula(paste("y ~ -1 + ", paste("x", 1:ncol(x), sep = ".", collapse = " + ")))
    dd <- data.frame(y = y , x = x)

    library("sgd")
    delta <- p / (2 * n)
    dd1 <- dd
    dd1$y <- (1 - 2 * delta) * dd1$y + delta
    fit_sgd <- sgd(form, data = dd1, model = "glm", model.control = list(family = "binomial"),
                   sgd.control = list(shuffle = TRUE, npasses = 20, method = "ai-sgd", lr = "adagrad"))


    dd <- simulate_candessur2020(n = 200000, kappa = 0.01, beta0 = 5, gamma = 20)
    ddf <- data.frame(Y = dd$Y, X = dd$X)
    form <- formula(paste("Y ~ -1 + ", paste("X", 1:ncol(dd$X), sep = ".", collapse = " + ")))



    system.time(
        m0 <- JeffreysMPL(y = y, m = NULL, X = x, a = 1/2, link = "logit", epsilon = 1e-04)
    )

    system.time(
        m1 <- glm(form, data = dd, family = binomial(), method = "brglm_fit", epsilon = 1e-04)
    )

    system.time(
        m2ml <- bigglm(form, data = ddf, family = binomial(), type = "ML", maxit = 100, epsilon = 1e-03, implementation = "1pass", chunksize = 1000, verbose = TRUE)
    )

    system.time(
        m2 <- bigglm(form, data = dd, family = binomial(), type = "BRASE", maxit = 100, tolerance = 1e-03, implementation = "1pass", chunksize = 1000, verbose = 1, start = coef(fit_sgd))
    )


    par(mfrow = c(1, 2))
    plot(coef(m2), coef(m0)); abline(0, 1)
    plot(coef(m2), coef(m1)); abline(0, 1)

    par(mfrow = c(1, 2))
    plot(coef(m2ml))
    points(beta, type ="l", col = "red", lwd =2)
    plot(coef(m2))
    points(beta, type ="l", col = "red", lwd =2)


}
