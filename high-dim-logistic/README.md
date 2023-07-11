## Distributed as part of the supplementary material for the manuscript
## "Bounded-memory adjusted scores estimation in generalized linear models with large data sets"
##
## Author: Ioannis Kosmidis
## Date: 09 July 2023
## Licence: GPL 2 or greater
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is".
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

This directory provides the data and scripts to reproduce the case
study for modelling Diverted US flights in 2000 in Section 4 of the
manuscript, and in Sections S2-S5 of the Supplementary Material
document.


The contributed R packages needed are dplyr, ggplot2, ggpp,
patchwork. The results have been obtained using R 4.3.1 with dplyr
v1.1.2, ggplot2 v3.4.2, ggpp v0.5.2, patchwork v1.1.2.


The easiest way to reproduce the relevant outputs in the manuscript
and Section S2 and Section S3 of the Supplementary Material document
is going through the following steps for each combination of n and
each of the two beta^* settings (see Section 5.2 of the main text for
what n and beta^* are, and Section S2 and Section S3 of the
Supplementary Material document, for a description of setting a and
setting b, respectively). The steps below are for n = 1000 and setting
a for beta^*.

1. edit "1-pt-simulate.R" and change the path in the
   `devtools::load_all` call to the path to the port of the biglm R
   package (also distributed as part of the supplementary).


2. In a terminal run 

> Rscript --no-init-file -e 'base_path <- "."; nobs <- 1000; b_setting <- "a"; ncores <- 5; source(file.path(base_path, "1-pt-simulate.R"))'

  where `base_path` should be set to the directory of this README
  file, and `ncores` is the number of cores you would like to use for
  the computation, `nobs` is n, and `b_setting` is the beta^* setting
  to consider. As the script is running the results/ directory will be
  populated with the results for each kappa, gamma and rho combination
  we consider in the experiment.

3. In a terminal run

> Rscript --no-init-file -e 'base_path <- "."; nobs <- 1000; b_setting <- "a"; source(file.path(base_path, "2-pt-plot.R"))'

  where the variables are set as in step 2. As the script is running
  the figures/ directory will be populated with the figures shown in
  Sections S2 and S3 in the Supplementary Materials document.


Steps 2 and 3 should be repeated for all remaining combinations of n
(`nobs <- 1000`, `nobs <- 2000`, `nobs <- 3000`, and beta^* settings
(`b_setting <- "a"`, `b_setting <- "b"`).


Once all results have been computed, the runtime summaries in Section
S4 and Section S5 can be reproduced by running in a terminal

> Rscript --no-init-file -e 'base_path <- "."; nobs <- 1000; b_setting <- "a"; rhosq <- 0; source(file.path(base_path, "3-pt-summarize-performance.R"))'
> Rscript --no-init-file -e 'base_path <- "."; nobs <- 1000; b_setting <- "b"; rhosq <- 0; source(file.path(base_path, "3-pt-summarize-performance.R"))'
> Rscript --no-init-file -e 'base_path <- "."; nobs <- 2000; b_setting <- "a"; rhosq <- 0; source(file.path(base_path, "3-pt-summarize-performance.R"))'
> Rscript --no-init-file -e 'base_path <- "."; nobs <- 2000; b_setting <- "b"; rhosq <- 0; source(file.path(base_path, "3-pt-summarize-performance.R"))'
> Rscript --no-init-file -e 'base_path <- "."; nobs <- 3000; b_setting <- "a"; rhosq <- 0; source(file.path(base_path, "3-pt-summarize-performance.R"))'
> Rscript --no-init-file -e 'base_path <- "."; nobs <- 3000; b_setting <- "b"; rhosq <- 0; source(file.path(base_path, "3-pt-summarize-performance.R"))'
