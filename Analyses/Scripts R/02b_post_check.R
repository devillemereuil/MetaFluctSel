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
options(mc.cores = ifelse(availableCores() > 4,
                          10,
                          availableCores()))

theme_set(theme_bw())

source("00_functions.R")

## ---------------- Performing the posterior predictive check

# Loading meta-data
mods <- readRDS("../STAN_mods.rds") %>%
        filter(Dist == "Poisson", Model_Year == "AR1")

# Performing the posterior predictive checks
restart_pb(nrow(mods))
mods <- 
    bind_cols(
        mods,
        pmap_dfr(list(mods[["ID"]],
                      mods[["Data"]]),
                 load_and_ppc,
                 progress    = TRUE,
                 plot        = FALSE)
    )

# Listing models needing zero-inflation
zi <-
    mods %>%
    filter(ppp_zero > 0.95) %>%
    chuck("ID")%>%
    str_remove("_AR1")

saveRDS(zi, file = "../need_zi.rds", version = 2)
