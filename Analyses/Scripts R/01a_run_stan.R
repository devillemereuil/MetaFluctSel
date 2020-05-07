library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(progress)
library(furrr)
plan(multiprocess)
filter  <- dplyr::filter
lag     <- dplyr::lag

theme_set(theme_bw())

source("00_functions.R")

nchains <- if_else(detectCores() > 10, 10, 4)
options(mc.cores = nchains)
rstan_options(auto_write = TRUE)

## ---------------- Loading the data

print(str_glue("Loading and formatting the data..."))

## Getting the data
alldata <- readRDS("../../Data/alldata_final.rds")
# NOTE Some df have quite a lot of missing IDs (up to 69%)

## Formatting the data for analysis (see 00_functions.R)
mods <- format_to_mods(alldata) %>%
        mutate(Dist = Dist %>%
                      recode(poisson    = "Poisson",
                             binomial   = "Binom"),
               Dist = if_else(map_lgl(Data, ~{min(.[["Fitness"]]) > 0}) & Dist == "Poisson",
                              "TruncPoisson",
                              Dist))

## ------------------ Setting up parameters for running the models

print(str_glue("Setting up everything..."))

## Setting up possible models for years
models_year <- c("FLAT", "EXPNULL", "EXPIID", "EXPAR1", "NULL", "IID", "AR1")
stanmods <-
    mods %>%
    select(Dist) %>%
    distinct() %>%
    crossing(tibble(Model_Year = models_year)) %T>%
    {restart_pb(nrow(.))} %>%
    mutate(Model = map(str_glue("../STAN models/model_dist{Dist}_type{Model_Year}.stan"),
                       ~ {pb$tick(0);pb$tick();stan_model(file = .)}))

## Now joining with the rest of mods
# This should prevent errors with loading too many dynamic libraries
mods <-
    mods %>%
    full_join(stanmods)

## ------------------ Running the models

print(str_glue("Fitting the models..."))

# Formatting some unique ID
mods <-
    mods %>%
    mutate(ID = str_glue("{DataID}_{Species}_{Population}_{Model_Year}"))

# Now running the models (fit_hmc saves the model to ../Output)
restart_pb(nrow(mods))
with(mods,
     pmap(list(Data, Model, ID),
          ~ fit_hmc(data    = ..1,
                    model   = ..2,
                    id      = ..3,
                    nchains = nchains,
                    progress= TRUE)))

# Keeping the output paths in the data.frame
mods <-
    mods %>%
    mutate(Model_Path = str_glue("../Output/{ID}_model.rds"))

## ------------------ Saving the tbl containing the model features

# Formatting the tbl to save it 
mods <-
    mods %>%
    select(DataID, Species, Population, ID, Data, Dist, Model_Year, Model_Path)

# Saving it
saveRDS(mods, file = "../STAN_mods.rds", version = 2)
