library(tidyverse)
library(magrittr)
filter  <- dplyr::filter
lag     <- dplyr::lag

source("00_functions.R")


## --------------------------- Getting the data
alldata <- readRDS("../../Data/alldata_summedFitness.rds")

## --------------------------- Formatting the data for analysis
mods <- format_to_mods(alldata)

## --------------------------- Now creating the IDs
ids <-
    mods %>%
    select(DataID, Species, Population) %>%
    distinct() %>%
    mutate(ID = str_replace(Species,
                            "([A-Z])[a-z]+ ([a-z]{2})[a-z ]+",
                            "\\1\\2")) %>%
    group_by(Species) %>%
    mutate(ID = if_else(rep(n() > 1, n()),
                        str_c(ID, 1:n()),
                        ID)) %>%
    ungroup() %>%
    mutate(ID = factor(ID,
                       levels = str_sort(ID, numeric = TRUE))) %>%
    select(ID, everything())

## ----------------------- Checking the order is the same as the summary table already generated
tbl_sum <- 
    read_csv("../../Data/Summary.csv") %>%
    select(ID, Species, Population)

# If the antijoin has no row, it means perfect identity
if (anti_join(tbl_sum, select(ids, -DataID)) %>% nrow() %>% equals(0)) {
    saveRDS(ids, "../FinalIDs.rds", version = 2)
} else {
    # Make sure we don't use wrong IDs: safer to delete them until this is sorted out!
    system("rm ../FinalIDs.rds")
    stop("IDs are not identical with the IDs in Data/Summary.csv, please check")
}


