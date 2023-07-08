# Simulate a the Candes and Sur (2020) setting
simulate_candessur2020 <- function(n, kappa = 0.2, gamma = 1, beta0 = 0, beta_star = rnorm(p)) {
    ## Section 2.4
    p <- ceiling(kappa * n)
    gamma0 <- sqrt(gamma^2 - beta0^2)
    beta <- gamma0 * beta_star / sqrt(sum(beta_star^2))
    X <- cbind(1, matrix(rnorm(n * p), n, p))
    eta <- drop(X %*% c(beta0, beta))
    Y <- 2 * (plogis(eta) > runif(n)) - 1
    out <- list(Y = Y, X = X)
    attr(out, "gamma0") <- gamma0
    attr(out, "beta") <- beta
    attr(out, "beta0") <- beta0
    attr(out, "kappa") <- kappa
    out
}

## Phase transition
compute_PT <- function(beta0 = 0,
                       gamma_grid = seq(0, 20, length = 100),
                       nsimu = 1000000,
                       ncores = 10,
                       XZU = NULL) {
    if (is.null(XZU)) {
        X <- rnorm(nsimu)
        U <- runif(nsimu)
        Z <- rnorm(nsimu)
    } else {
        X <- XZU$X
        U <- XZU$U
        Z <- XZU$Z
    }
    kappa <- function(gamma) {
        gamma0 <- sqrt(gamma^2 - beta0^2)
        Y <- -1 + 2 * (plogis(beta0 + gamma0 * X) > U)
        V <- X * Y
        obj <- function(ts) {
            mean(pmax(ts[1] * Y + ts[2] * V - Z, 0)^2)
        }
        optim(c(0, 0), obj, method = "BFGS")$value
    }
    kappa <- unlist(mclapply(1:length(gamma_grid), function(k) {
        kappa(gamma = gamma_grid[k])
    }, mc.cores = ncores))
    data.frame(kappa = kappa, gamma = gamma_grid)
}

get_results <- function(setting,
                        beta0 = 0,
                        nz_perc = 0.2,
                        beta_star_setting = "a",
                        verbose = 1,
                        ncores = 2) {
    kappa <- setting$kappa
    gamma <- setting$gamma
    seed <- setting$seed
    maxit <- setting$maxit
    tolerance <- setting$tolerance
    chunksize <- setting$chunksize

    n <- setting$n
    p <- setting$p
    nz <- ceiling(p * nz_perc) ## nz_perc -x and nz_perc +x
    beta_star <- switch(beta_star_setting,
                        "a" = seq(-10, 10, length.out = p),
                        "b" = c(rep(-10, nz), rep(10, nz), rep(0, p - 2 * nz)),
                        stop("invalid beta star setting"))
    set.seed(seed)

    datasets <- lapply(1:setting$repetitions, function(j) {
        dat <- simulate_candessur2020(n = n,
                                      kappa = kappa,
                                      gamma = gamma,
                                      beta0 = beta0,
                                      beta_star = beta_star)
        dd <- data.frame(Y = dat$Y, X = dat$X)
        dd$Y <- 0.5 * (dd$Y + 1)
        attr(dd, "beta0") <- attr(dat, "beta0")
        attr(dd, "beta") <- attr(dat, "beta")
        dd
    })

    form <- as.formula(paste("Y ~ -1 + ", paste("X", 1:(p + 1), sep = ".", collapse = " + ")))

    results <- mclapply(1:setting$repetitions, function(j) {
        dd <- datasets[[j]]
        true_betas <- c(attr(dd, "beta0"), attr(dd, "beta"))
        time_br <- system.time(
            fit_br <- bigglm(form, family = binomial(), data = dd,
                             type = "BRASE",
                             implementation = "2pass",
                             tolerance = tolerance, maxit = maxit,
                             verbose = verbose,
                             chunksize = chunksize)
        )
        coefs_br <- coef(fit_br)
        ses_br <- summary(fit_br)$mat[, "SE"]
        elapsed_br <- time_br[["elapsed"]]
        iter_br <- fit_br$iter
        libeta_br <- max(abs(fit_br$delta))

        if (setting$mle_exists) {
            time_ml <- system.time(
                fit_ml <- bigglm(form, family = binomial(), data = dd,
                                 verbose = verbose,
                                 tolerance = tolerance,
                                 maxit = maxit,
                                 chunksize = chunksize,
                                 start = coefs_br)
            )
            elapsed_ml <- time_ml[["elapsed"]]
            iter_ml <- fit_ml$iter
            libeta_ml <- max(abs(fit_ml$delta))
            ## If ML did not converge then set estimates and ses to NA
            if (libeta_ml > 10) {
                coefs_ml <- ses_ml <- rep(NA, length(true_betas))
            } else {
                coefs_ml <- coef(fit_ml)
                ses_ml <- summary(fit_ml)$mat[, "SE"]
            }
        } else {
            coefs_ml <- ses_ml <- rep(NA, length(true_betas))
            elapsed_ml <- iter_ml <- libeta_ml <- NA
        }


        cat("Setting:", "n =", n,
            "| κ =", kappa, "| γ =", gamma,
            "rep: ", j, "/", setting$repetitions,
            "| |delta_ML|_inf:", round(libeta_ml, 4), paste0("(", iter_ml,")"),
            "| |delta_BR|_inf:", round(max(libeta_br), 4), paste0("(", iter_br,")"),"\n")

        results <- data.frame(method = rep(c("ML", "mJPL"), each = length(true_betas)),
                              estimate = c(coefs_ml, coefs_br),
                              se = c(ses_ml, ses_br),
                              truth = c(true_betas, true_betas),
                              parameter = rep(seq_along(true_betas) - 1, 2),
                              kappa = kappa,
                              gamma = gamma,
                              rhosq = setting$rhosq,
                              repetition = j,
                              mle_exists = setting$mle_exists,
                              boundary = setting$boundary,
                              n = n,
                              p = p)
        performance <- data.frame(method = c("ML", "mJPL"),
                                  elapsed = c(elapsed_ml, elapsed_br),
                                  iter = c(iter_ml, iter_br),
                                  libeta = c(libeta_ml, libeta_br),
                                  kappa = kappa,
                                  gamma = gamma,
                                  rhosq = setting$rhosq,
                                  repetition = j,
                                  mle_exits = setting$mle_exists,
                                  boundary = setting$boundary,
                                  n = n,
                                  p = p)
        rownames(results) <- 1:nrow(results)
        list(results = results, performance = performance)
    }, mc.cores = ncores)

    performance <- do.call("rbind", lapply(results, "[[", "performance"))
    results <- do.call("rbind", lapply(results, "[[", "results"))
    attr(results, "performance") <- performance
    attr(results, "formula") <- form
    attr(results, "seed") <- seed
    results
}

correct_mJPL_estimates <- function(results) {
    results |> mutate(estimate = ifelse(!boundary & method == "mJPL" & !mle_exists,
                                        estimate * kappa * gamma / (1 - rhosq)^(1/2),
                                        estimate))
}

summarize_results <- function(results, across_parameters = FALSE) {
    if (across_parameters) {
        results <- results |>
            group_by(method, kappa, gamma, rhosq, mle_exists)
    } else {
        results <- results |>
            group_by(method, parameter, kappa, gamma, rhosq, mle_exists)
    }
    results |>
        summarize(estimated_mean = mean(estimate),
                  bias = mean(estimate - truth),
                  rmse = sqrt(mean((estimate - truth)^2)),
                  mad = mean(abs(estimate - truth)),
                  pu = mean(estimate < truth),
                  repetitions = n(),
                  estimated_sd = sd(estimate),
                  avg_se = mean(se)) |>
        data.frame()

}

## For a single kappa/gamma combination
plot_results <- function(summaries, cols = c("#CF4446", "#00AD9A", "#FB9A06"),
                         p_alpha = 0.2, p_size = 1,
                         type = c("estimate_vs_truth", "estimate_and_truth")) {
    type <- match.arg(type)
    summaries$method <- factor(summaries$method, levels = c("ML", "mJPL"), ordered = TRUE)
    ## If the ML estimates do not exist plot nothing
    if (all(!summaries$mle_exists)) {
        summaries <- summaries |> subset(method != "ML")

    }
    lims <- with(summaries, range(c(truth, estimate)))
    if (type == "estimate_vs_truth") {
        p1 <- ggplot(summaries) +
            geom_point(aes(x = truth, y = estimate, col = method), alpha = p_alpha, size = p_size) +
            geom_abline(aes(intercept = 0, slope = 1), col = "black", lwd = 0.5) +
            geom_smooth(aes(x = truth, y = estimate),
                        method = "lm", formula = "y ~ x",
                        se = FALSE, col = cols[2], lwd = 0.5) +
            coord_cartesian(x = lims, y = lims) +
            scale_colour_manual(values = c("ML" = cols[1], "mJPL" = cols[3]))
    }
    if (type == "estimate_and_truth") {
        sd <- summaries |> group_by(kappa, gamma, method, truth) |>
            summarize(e_mean = mean(estimate)) |> inner_join(summaries, c("kappa", "gamma", "method", "truth"))
        p1 <- ggplot(summaries) +
            geom_point(aes(x = parameter, y = estimate, col = method), alpha = p_alpha, size = p_size) +
            geom_line(aes(x = parameter, y = truth), col = "black", lwd = 0.5) +
            geom_line(data = sd, aes(x = parameter, y = e_mean), col = cols[2], lwd = 0.5) +
           coord_cartesian(y = lims) +
            scale_colour_manual(values = c("ML" = cols[1], "mJPL" = cols[3]))
    }
    p1 + theme_minimal() +
        facet_grid(~ method) +
        theme(legend.position = "none") +
        theme(
            plot.background = element_rect(fill= 'transparent', color = "grey"),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_blank()
        ) +
        labs(x = NULL, y = NULL)
}



plot_PT <- function(PT, max_kappa = NULL, cols = c("#F1F1F1", "#FFFFFF")) {
    ## cols = c("#ffd1e4", "#5ff9ff")) {
    PT <- PT[order(PT$kappa), ]
    max_gamma <- max(PT$gamma)
    max_kappa <- ifelse(is.null(max_kappa), max(PT$kappa), max_kappa)
    polygon1 <- data.frame(kappa = c(PT$kappa, max_kappa, max_kappa),
                           gamma = c(PT$gamma, 0, max_gamma))
    polygon2 <- data.frame(kappa = c(0, 0, PT$kappa),
                           gamma = c(0, max_gamma, PT$gamma))
    ggplot(PT) +
        geom_polygon(data = polygon1, aes(kappa, gamma), fill = cols[1]) +
        geom_polygon(data = polygon2, aes(kappa, gamma), fill = cols[2]) +
        geom_line(aes(kappa, gamma), col = "grey") +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(x = expression(kappa), y = expression(gamma))
}




