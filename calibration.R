### Script to calibrate the EAAA P.a. model to the target features as listed in the MaxART modelling report.

# 0. Loading libraries
library(RSimpactCyan)
library(devtools)
install_github("wdelva/RSimpactHelp")
library(RSimpactHelper)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(metafolio)
library(EasyABC)


# 1. Target features
features.pop.growth <- exp(0.015)
features.hiv.prev <- c(0.143, 0.008, 0.315, 0.066, 0.467, 0.213, 0.538, 0.366, 0.491, 0.47, 0.397, 0.455, 0.316, 0.425)
features.hiv.inc <- exp(c(0.038, 0.008, 0.043, 0.016, 0.02, 0.026, 0.027, 0.031, 0.04, 0.004, 0.021, 0.012, 0.012, 0))
features.art.cov <- c(0.37, 0.40, 0.44, 0.49, 0.58, 0.67, 0.76, 0.85) # c(0.33, 0.38, 0.45, 0.51, 0.61, 0.7, 0.8)
features.vl.suppr <- 0.74  # 0.68
target.features.EAAA <- c(features.pop.growth, features.hiv.prev, features.hiv.inc, features.art.cov, features.vl.suppr)


# 2. Prior distributions for parameters to be calibrated
priors.EAAA <-list(c("unif", 1.0, 1.5), # hivtransmission.param.f1 = 1.1
                   c("unif", 0.0, 0.5), # formation.hazard.agegapry.gap_agescale_man and ~_woman = 0.25
                   c("unif", -3.0, 3.0), # person.agegap.man.dist.normal.mu and ~.woman.~ = 0
                   c("unif", 1.0, 5.0), # person.agegap.man.dist.normal.sigma and ~.woman.~ = 2.5
                   c("unif", 0, 1.5), # person.eagerness.man.dist.gamma.a = 0.23 # Used to be 0.1, 1.5
                   c("unif", 0, 1.5), # person.eagerness.woman.dist.gamma.a = 0.23 # Used to be 0.1, 1.5
                   c("unif", 0, 80), # person.eagerness.man.dist.gamma.b = 45 # Used to be 10, 80
                   c("unif", 0, 80), # person.eagerness.woman.dist.gamma.b = 45 # Used to be 10, 80
                   c("unif", -2, -0.2), # formation.hazard.agegapry.gap_factor_man_exp and ~_woman_~ = -0.7
                   c("unif", -1.0, 5.0), # formation.hazard.agegapry.baseline = 2.8
                   c("unif", -1, -0.05), # formation.hazard.agegapry.numrel_man = -0.5
                   c("unif", -1, -0.05), # formation.hazard.agegapry.numrel_woman = -0.5
                   c("unif", -5, -1.0), # conception.alpha_base = -2.7
                   c("unif", -2, 1), # dissolution.alpha_0 = -0.52
                   c("unif", -7, -3), #  art.intro["diagnosis.baseline"] = -2
                   c("unif", 0, 4), # jump to art.intro1["diagnosis.baseline"] # previously 0.2
                   c("unif", 0, 4), # jump to art.intro1["diagnosis.baseline"] # previously 0.3
                   c("unif", 0, 4), # jump to art.intro1["diagnosis.baseline"] # previously 0.5
                   c("unif", -2, 2)) # jump to art.intro1["diagnosis.baseline"] # previously 0



inputvector.EAAA.example <- c(12345,
                              1.1,
                              0.4,
                              1,
                              3.5,
                              0.8,
                              0.3,
                              45,
                              30,
                              -1,
                              3.8,
                              -0.5,
                              -0.45,
                              -2.5,
                              -1,
                              -5,
                              2,
                              2,
                              1,
                              0)

testing.EAAA.calibration <- EAAA.calibration.wrapper(inputvector = inputvector.EAAA.example)


prior.boundaries.booleans <- !unlist(priors.EAAA) %in% "unif"
boundaries.matrix <- unlist(priors.EAAA)[prior.boundaries.booleans] %>%
  as.numeric() %>% matrix(byrow = TRUE,
                          ncol = 2)

# 3. Calibration with Sequential ABC (Lenormand's Adaptive Population Monte Carlo ABC)
EAAA.Seq <- ABC_sequential(model = EAAA.calibration.wrapper,
                           method = "Lenormand",
                           prior = priors.EAAA,
                           summary_stat_target = target.features.EAAA,
                           nb_simul = 1000,
                           alpha = 0.25,
                           p_acc_min = 0.5,
                           use_seed = TRUE,
                           seed_count = 0,
                           n_cluster = 4,
                           inside_prior = TRUE)

