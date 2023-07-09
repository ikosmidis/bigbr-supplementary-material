library("memisc")
devtools::load_all("~/Repositories/biglm")

if (interactive()) {
    experiment_path <- "~/Repositories/bigbr_supplementary_material/diverted-flights/"
    data_path <- file.path(experiment_path, "data")
    results_path <- file.path(experiment_path, "results")
    load(file.path(results_path, "diverted-fits.rda"))
}

## Helper functions to interface memisc::mtable
getSummary.biglm <- function(obj, alpha = 0.05) {
    x <- summary(obj)$mat
    coefs <- x[, "Coef"]
    ses <- x[, "SE"]
    stat <- coefs/ses
    p <- 2 * pnorm(abs(stat), lower.tail = FALSE)
    cc <- cbind(coefs, ses, stat, p)
    colnames(cc) <- c("est", "se", "stat", "p")
    rownames(cc) <- rownames(x)
    list(coef = cc, call = obj$call,
         time = structure(obj$time["elapsed"],
                          names = "Time (sec)"),
         iterations = structure(obj$iterations,
                                names = "Iterations"),
         time_per_iteration = structure(obj$time["elapsed"] / obj$iterations,
                                        names = "Time / Iteration (sec)"),
         deviance = structure(deviance(results[[1]]), names = "Deviance"))
}

## Set names to parameters as in expression (12) of the paper
covariate_names <- results[[1]]$names
parameter_names <- c("$\\alpha$",
                     paste0("$\\beta_", 2:12, "$"),
                     paste0("$\\gamma_", 2:7, "$"),
                     paste0("$\\delta_", 2:11, "$"),
                     "$\\zeta_{(d)}$",
                     "$\\zeta_{(a)}$",
                     "$\\rho$",
                     paste0("\\psi_{(d),", 1:3, "}$"),
                     paste0("\\psi_{(a),", 1:3, "}$"))
results <- lapply(results, function(x) {
    x$names <- parameter_names
    x
})

## Parameter name - Coveriate names correspondence
cbind(covariate_names, parameter_names)

memisc:::mtable("ML-1pass" = results[["ML-1pass-15"]],
                "ML-1pass" = results[["ML-1pass-20"]],
                "BRASE-1pass" = results[["BRASE-1pass-20"]],
                "BRASE-2pass" = results[["BRASE-2pass-20"]],
                "MJPL-1pass" = results[["MJPL-1pass-20"]],
                "MJPL-2pass" = results[["MJPL-2pass-20"]],
                digits = 2,
                signif.symbols = NULL,
                summary.stats = c("")) |> toLatex()

