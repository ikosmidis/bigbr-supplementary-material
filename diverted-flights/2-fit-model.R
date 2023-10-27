if (interactive()) {
    devtools::load_all("~/Repositories/biglm")
    library("parallel")
    experiment_path <- "~/Repositories/bigbr-supplementary-material/diverted-flights/"
    data_path <- file.path(experiment_path, "data")
    results_path <- file.path(experiment_path, "results")
    air <- readRDS(file.path(data_path, "air2000_combined.rds"))
    n_cores <- 1
}

form_air <- Diverted ~ Month +  DayOfWeek + UniqueCarrier + CRSDepTime + CRSArrTime + Distance + Orig_x + Orig_y + Orig_z + Dest_x + Dest_y + Dest_z

## family object
fam <- binomial("probit")

run_settings <- data.frame(maxit = c(15, 20, 20, 20, 20, 20),
                           tolerance = 1e-03,
                           chunksize = 1e+05,
                           type = c("ML", "ML", "BRASE", "BRASE", "MJPL", "MJPL"),
                           implementation = c("1pass", "1pass", "1pass", "2pass", "1pass", "2pass"),
                           verbose = FALSE)
row.names(run_settings) <- with(run_settings, paste(type, implementation, maxit, sep = "-"))

results <- mclapply(1:nrow(run_settings), function(s) {
    timing <- system.time(
        mod <- with(run_settings[s, ],
                    bigglm(form_air, data = air, family = fam,
                           maxit = maxit,
                           tolerance = tolerance,
                           chunksize = chunksize,
                           type = type,
                           implementation = implementation,
                           verbose = verbose)
                    )
    )
    cat(mod$type, mod$implementation, "complete in")
    print(timing)
    mod$time <- timing
    mod
}, mc.cores = n_cores)
names(results) <- row.names(run_settings)

## mBR and mJPL fits using brglm2
library("brglm2")
timing_mBR <- system.time(
    mod_mBR <- glm(form_air, data = air, family = fam,
                   maxit = 20,
                   tolerance = 1e-03,
                   method = "brglm_fit",
                   type = "AS_mean",
                   max_step_factor = 1,
                   start = rep(0, 37))
)

timing_mJPL <- system.time(
    mod_mJPL <- glm(form_air, data = air, family = fam,
                   maxit = 20,
                   tolerance = 1e-03,
                   method = "brglm_fit",
                   type = "MPL_Jeffreys",
                   max_step_factor = 1,
                   start = rep(0, 37))
)


save(results, run_settings, form_air, timing_mBR, timing_mJPL,
     file = file.path(results_path, "diverted-fits.rda"))

