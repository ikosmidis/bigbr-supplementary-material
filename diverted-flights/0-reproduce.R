devtools::load_all("~/Repositories/biglm")
library("parallel")
library("memisc")
library("dplyr")

n_cores <- 1
experiment_path <- "~/Repositories/bigbr-supplementary-material/diverted-flights/"
data_path <- file.path(experiment_path, "data")
results_path <- file.path(experiment_path, "results")

## Prepare data
cat("Preparing data ...")
air <- read.csv(file.path(data_path, "2000.csv.bz2"))
airports <- read.csv(file.path(data_path, "airports.csv"))
source(file.path(experiment_path, "1-prepare-data.R"))
cat("Done!\n")

## Fit model
cat("Fitting the models ...")
source(file.path(experiment_path, "2-fit-model.R"))
cat("Done!\n")

## Outputs
cat("Returning the outputs ...")
source(file.path(experiment_path, "3-output-tables.R"))
cat("Done!\n")
