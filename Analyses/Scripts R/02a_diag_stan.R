library(tidyverse)
library(magrittr)
library(rstan)
library(coda)
library(bayesplot)
library(furrr)
library(progress)
filter  <- dplyr::filter
lag     <- dplyr::lag
extract <- rstan::extract

plan(multiprocess)
options(mc.cores = ifelse(availableCores() > 10,
                          10,
                          availableCores() - 1))

theme_set(theme_bw())

source("00_functions.R")

## ---------------- Performing the diagnostics

# Loading meta-data
mods <- readRDS("../STAN_mods.rds")

# Performing the diagnostics
mods <- 
    bind_cols(
        mods,
        future_pmap_dfr(list(mods[["ID"]],
                             mods[["Model_Year"]],
                             mods[["Dist"]]),
                        load_and_diagnose,
                        progress    = FALSE,
                        trace       = FALSE,        # NOTE Set to TRUE for first run
                        .progress   = TRUE)
    )


## ---------------- Checking some diagnostics

mods <-
    mods %>%
    mutate(Issue_STAN = (map(Diag_STAN, "Divergence") %>% {. > 100} |
                         map(Diag_STAN, "Energy") %>% {. > 0}))

mods <-
    mods %>%
    mutate(Issue_MCMC = (map(Diag_MCMC, "N_eff_ratio") %>% map_dbl(min) %>% {. < 0.05}) |
                        (map(Diag_MCMC, "Rhat") %>% map_dbl(max) %>% {. > 1.05}))
