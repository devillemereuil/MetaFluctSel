library(tidyverse)
library(magrittr)
library(rstan)
library(loo)
library(progress)
library(furrr)
plan(multisession)
options(mc.cores = ifelse(availableCores() > 10,
                          10,
                          availableCores() - 1))
filter  <- dplyr::filter
lag     <- dplyr::lag
extract <- rstan::extract

theme_set(theme_bw())

source("00_functions.R")

## ---------------- Loading the meta-data

print(str_glue("Loading the meta-data"))

# Loading meta-data
mods <- readRDS("../STAN_mods.rds") %>%
        # NOTE Removing the NULL models, replaced by VARINT
        filter(!str_detect(Model_Year, "NULL"))
ids  <- readRDS("../FinalIDs.rds")

## ---------------- Computing Likelihood and performing LOO

print(str_glue("Performing the LOO computation"))

restart_pb(nrow(mods))
mods <-
    mods %>%
    mutate(LOO = future_map2(ID, Data, load_and_loo,
                             progress = FALSE,
                             .progress = TRUE))             # .progress is a future_* argument

## ---------------- Extracting LOOIC and p_loo

print(str_glue("Extracting information and formating the output"))
mods <-
    mods %>%
    mutate(LOOIC = map(LOO, "estimates") %>% map_dbl(~ .["looic", "Estimate"]),
           P_loo = map(LOO, "estimates") %>% map_dbl(~ .["p_loo", "Estimate"]))

## Extracting and formatting LOOIC information
loo <-
    mods %>%
    select(DataID, Species, Population, Dist, Model_Year, LOOIC) %>%
    group_by(DataID, Species, Population) %>%
    mutate(Delta_LOOIC = LOOIC - min(LOOIC),
           LOO_Weight = exp(-0.5 * Delta_LOOIC) / sum(exp(-0.5 * Delta_LOOIC))) %>%
    arrange(Delta_LOOIC, .by_group = TRUE) %>%
    ungroup()

print(str_glue("Saving and done!"))
saveRDS(loo, file = "../loo.rds", version = 2)

## Making a "square" version of LOOIC analysis with 4 models
loo_square <-
    loo %>%
    filter(str_detect(Model_Year, "AR1|VARINT")) %>%
    select(DataID, Species, Population, Dist, Model_Year, LOOIC) %>%
    group_by(DataID, Species, Population) %>%
    mutate(Delta_LOOIC = LOOIC - min(LOOIC),
           LOO_Weight = exp(-0.5 * Delta_LOOIC) / sum(exp(-0.5 * Delta_LOOIC))) %>%
    arrange(Delta_LOOIC, .by_group = TRUE) %>%
    ungroup()

## ---------------- Creating the table for the article
out <-
    loo %>%
    left_join(ids) %>%
    arrange(ID) %>%
    select(ID, everything(), -DataID, -Population, -Species) %>%
    group_by(ID) %>%
    ungroup()
write_csv(out, path = "../loo.csv")

## ---------------- Summarising the results

## Loading information about taxa
taxa <- readRDS("../../Data/taxa.rds")

## Best Models
best <-
    loo %>%
    filter(Delta_LOOIC == 0) %>%
    left_join(taxa)

sum_best <-
    best %>%
    group_by(Taxon) %>%
    summarise(Table = list(table(Model_Year))) %>%
    mutate(Table = map(Table, enframe)) %>%
    unnest_legacy() %>%
    complete(name, nesting(Taxon), fill = list(value = 0)) %>%
    spread(name, value)

write_csv(sum_best, path = "../loo_best_taxa.csv")

## Relative support
rel <-
    loo %>%
    # Adding the taxon information
    left_join(taxa) %>%
    # Adding the whole dataset against, but Taxon = total
    bind_rows(loo %>% mutate(Taxon = "Total")) %>%
    # Group by the taxa
    group_by(Taxon) %>%
    # Compute total weight
    mutate(Tot_Weight = sum(LOO_Weight)) %>%
    # Group by models and taxa
    group_by(Model_Year, Taxon) %>%
    # Compute relative support
    summarise(Rel_Sup = sum(LOO_Weight)/unique(Tot_Weight)) %>%
    # Spread for a better format for article display
    spread(Model_Year, Rel_Sup)

write_csv(rel, path = "../loo_relsupport.csv")

## Relative support (article format)

generate_rel_table <-
    . %>%
    mutate(ID = str_replace(Species,
                            "([A-Z])[a-z]+ ([a-z]{2})[a-z ]+",
                            "\\1\\2")) %>%
    group_by(Species) %>%
    mutate(ID = if_else(rep(n() > 1, n()),
                        str_c(ID, 1:n()),
                        ID)) %>%
    ungroup() %>%
    left_join(taxa) %>%
    bind_rows(mutate(., Taxon = "Total")) %>%
    group_by(Taxon) %>%
    mutate(Tot_Weight = sum(LOO_Weight)) %>%
    group_by(Model_Year, Taxon) %>%
    summarise(Rel_Sup = sum(LOO_Weight)/unique(Tot_Weight)) %>%
    ungroup() %>%
    mutate(Model_Year = recode(Model_Year,
                               FLAT     = "NoSel",
                               VARINT   = "ConstOpt",
                               IID      = "FluctOpt",
                               AR1      = "FluctCorrOpt",
                               EXPVARINT= "ConstDir",
                               EXPIID   = "FluctDir",
                               EXPAR1   = "FluctCorrDir"),
           Model_Year = factor(Model_Year,
                               levels = c("NoSel",
                                          "ConstDir", "FluctDir", "FluctCorrDir",
                                          "ConstOpt", "FluctOpt", "FluctCorrOpt"))) %>%
    spread(Model_Year, Rel_Sup)

rel <- generate_rel_table(loo)
write_csv(rel, path = "../loo_relsupport_article.csv")

rel_square <- generate_rel_table(loo_square)
write_csv(rel_square, path = "../loo_relsupport_square.csv")

## Bigger table with models description

desc <-
    tribble(
        ~ ID,           ~ Shape,        ~ `Fluctuations`,    ~ `Autocorrelation`,
        "NoSel",        "Flat",         "\\ding{56}",        "\\ding{56}",
        "ConstDir",     "Monotonic",    "\\ding{56}",        "\\ding{56}",
        "ConstOpt",     "Gaussian",     "\\ding{56}",        "\\ding{56}",
        "FluctDir",     "Monotonic",    "\\ding{52}",        "\\ding{56}",
        "FluctOpt",     "Gaussian",     "\\ding{52}",        "\\ding{56}",
        "FluctCorrDir", "Monotonic",    "\\ding{52}",        "\\ding{52}",
        "FluctCorrOpt", "Gaussian",     "\\ding{52}",        "\\ding{52}"
    )

desc <-
    desc %>%
    left_join(rel %>%
              pivot_longer(-Taxon, names_to = "ID", values_to = "Weight") %>%
              pivot_wider(names_from = "Taxon", values_from = "Weight"))
write_csv(desc, path = "../loo_rel_with_description.csv")


## Create a graph for a talk

graph_rel <-
    rel %>%
    gather("Model", "Support", -Taxon) %>%
    mutate(Model = recode(Model,
                          FLAT     = "NoSel",
                          VARINT   = "ConstOpt",
                          IID      = "FluctOpt",
                          AR1      = "FluctCorrOpt",
                          EXPVARINT= "ConstDir",
                          EXPIID   = "FluctDir",
                          EXPAR1   = "FluctCorrDir"),
           Model = factor(Model,
                          levels = c("NoSel",
                                     "ConstDir", "FluctDir", "FluctCorrDir",
                                     "ConstOpt", "FluctOpt", "FluctCorrOpt")))

p <- ggplot(graph_rel) +
     geom_col(aes(x = Taxon, y = Support, fill = Model), colour = "black", position = "stack") +
     theme(text = element_text(family = "Noto Sans", size = 16))
cairo_pdf("../../Figures/loo_rel.pdf")
p
dev.off()
