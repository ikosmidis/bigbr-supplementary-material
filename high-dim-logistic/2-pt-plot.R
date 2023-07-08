## Rscript --no-init-file -e 'nobs <- 500; beta_star_setting <- "a"; source("2-pt-plot.R")'
## Rscript --no-init-file -e 'nobs <- 500; beta_star_setting <- "b"; source("2-pt-plot.R")'
## Rscript --no-init-file -e 'nobs <- 2000; beta_star_setting <- "a"; source("2-pt-plot.R")'
## Rscript --no-init-file -e 'nobs <- 2000; beta_star_setting <- "b"; source("2-pt-plot.R")'

devtools::load_all("~/Repositories/biglm")
library("ggplot2")
library("ggpp")
library("patchwork")
library("dplyr")
library("tibble")

options(dplyr.summarise.inform = FALSE)

## beta_star_setting and nobs are set by the following lines only if R
## is used interactively
if (interactive()) {
    nobs <- 2000
    beta_star_setting <- "a"
}

include_title <- TRUE
plot_type <- if (beta_star_setting == "a") "estimate_vs_truth" else "estimate_and_truth"
rhosq_grid <- c(0.0, 0.25, 0.5, 0.75, 0.8, 0.9)

exp_h <- 5 * 200 * 1.5
vp_h <- 0.145
exp_w <- exp_h * sqrt(2)
exp_ratio <- exp_h / exp_w
vp_w <- vp_h * exp_ratio
## transparency of the points
p_size <- 0.5
alpha_fac <- 150

code_path <- "~/Repositories/bigbr-supplementary-material/high-dim-logistic/"
results_path <- file.path(code_path, "results")
figures_path <- file.path(code_path, "figures")
source(file.path(code_path, "functions.R"))

for (rhosq in rhosq_grid) {

    base_name <- paste0("mJPL-setting-", beta_star_setting,
                        "-n-", nobs, "-rho2-",
                        format(rhosq, digits = 2, nsmall = 2))

    ## Get phase transfition curve for rhosq
    load(file.path(results_path, paste0("PT-n-", nobs, "-setting-", beta_star_setting, "-rhosq-", rhosq, ".rda")))
    ## Get estimates for rhosq
    files <- dir(results_path,
                 pattern = paste0("estimates-n-", nobs, "-setting-", beta_star_setting, "-rhosq-", rhosq, "-"),
                 full.names = TRUE)
    res <- NULL
    for (f in files) {
        load(f)
        res <- rbind(res, results)
    }

    ## Avoid any unexpectedly large or NA ML estimates (due to
    ## non-existence close to the PT curve)
    large_ML <- res$method == "ML" & abs(res$se) > 100
    NA_ML <- res$method == "ML" & is.na(res$estimate)
    res <- res[!(large_ML | NA_ML), ]

    ## Find unique kappa-gamma in estimates
    kappa_gamma <- unique(res[c("kappa", "gamma", "mle_exists")])

    ## Plot mJPL and ML estimates (the latter only when they exist)
    insets_estimates <- as.list(rep(NA, nrow(kappa_gamma)))
    for (wh in 1:nrow(kappa_gamma)) {
        ckappa <- kappa_gamma[wh, "kappa"]
        cgamma <- kappa_gamma[wh, "gamma"]
        p_alpha <- alpha_fac * (1 - ckappa) / nobs
        ep <- res |>
            dplyr::filter(parameter > 1, kappa == ckappa, gamma == cgamma) |>
            plot_results(p_alpha = p_alpha, p_size = p_size, type = plot_type)
        insets_estimates[[wh]] <- tibble(x = ckappa,
                                         y = cgamma,
                                         plot = list(ep))
    }
    out <- plot_PT(pt, max_kappa = 0.6) +
        geom_point(data = kappa_gamma, aes(x = kappa, y = gamma), pch = 18, size = 1,
                   col = "grey")

    for (wh in 1:nrow(kappa_gamma)) {
        ckappa <- kappa_gamma[wh, "kappa"]
        cgamma <- kappa_gamma[wh, "gamma"]
        out <- out +  geom_plot(data = insets_estimates[[wh]],
                                aes(x = x, y = y, label = plot),
                                vp.width = ifelse(kappa_gamma[wh, "mle_exists"], 2 * vp_w, vp_w),
                                vp.height = vp_h,
                                vjust = "bottom", hjust = "left")
    }
    tt <- substitute(paste(n == x0, ", ", rho^2 == y0),
                     list(x0 = nobs, y0 = rhosq))
    if (include_title) {
        out_u <- out + labs(title = tt)
    }
    ## png(file.path(figures_path, paste0(base_name, "_uncorrected.png")),
    ##     width = exp_w, height = exp_h, res = 300)
    ## print(out_u)
    ## dev.off()

    ## Plot mJPL and ML estimates (the latter only when they exist),
    ## adjusting the former by kappa * gamma /sqrt(1 - rhosq) when the
    ## ML estimates do not exist
    insets_estimates <- as.list(rep(NA, nrow(kappa_gamma)))
    for (wh in 1:nrow(kappa_gamma)) {
        ckappa <- kappa_gamma[wh, "kappa"]
        cgamma <- kappa_gamma[wh, "gamma"]
        p_alpha <- alpha_fac * (1 - ckappa) / nobs
        ep <- res |>
            correct_mJPL_estimates() |>
            dplyr::filter(parameter > 1, kappa == ckappa, gamma == cgamma) |>
            plot_results(p_alpha = p_alpha, p_size = p_size, type = plot_type)
        insets_estimates[[wh]] <- tibble(x = ckappa,
                                         y = cgamma,
                                         plot = list(ep))
    }
    out <- plot_PT(pt, max_kappa = 0.6) +
        geom_point(data = kappa_gamma, aes(x = kappa, y = gamma), pch = 18, size = 1,
                   col = "grey")
    for (wh in 1:nrow(kappa_gamma)) {
        ckappa <- kappa_gamma[wh, "kappa"]
        cgamma <- kappa_gamma[wh, "gamma"]
        out <- out +  geom_plot(data = insets_estimates[[wh]],
                                aes(x = x, y = y, label = plot),
                                vp.width = ifelse(kappa_gamma[wh, "mle_exists"], 2 * vp_w, vp_w),
                                vp.height = vp_h,
                                vjust = "bottom", hjust = "left")
    }
    tt <- substitute(paste(n == x0, ", ", rho^2 == y0),
                     list(x0 = nobs, y0 = rhosq))
    if (include_title) {
        out_c <- out + labs(title = tt)
    }
    ## png(file.path(figures_path, paste0(base_name, "_corrected.png")),
    ##     width = exp_w, height = exp_h, res = 300)
    ## print(out_c)
    ## dev.off()

    f_u <- out_u +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(x = NULL)
    f_c <- out_c + labs(title = NULL)

    png(file.path(figures_path, paste0(base_name, ".png")),
        width = exp_w, height = exp_h * 1.8, res = 300)
    print(f_u / f_c)
    dev.off()

 }


