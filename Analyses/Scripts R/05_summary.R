library(tidyverse)
library(magrittr)
library(rstan)
library(progress)
library(furrr)
library(patchwork)
plan(multiprocess)
options(mc.cores = ifelse(availableCores() > 8,
                          10,
                          availableCores()))
filter  <- dplyr::filter
lag     <- dplyr::lag
extract <- rstan::extract

theme_set(theme_bw())

source("00_functions.R")

## ---------------- Loading the models

# Loading meta-data
mods <- readRDS("../STAN_mods.rds") 
loo  <- readRDS("../loo.rds")
ids  <- readRDS("../FinalIDs.rds")

# Subsetting the two most complete models
mods <-
    mods %>%
    left_join(loo %>% select(DataID, Species, Population, Model_Year, Delta_LOOIC, LOO_Weight)) %>%
    group_by(DataID, Species, Population) %>%
    filter(!any(Model_Year == "FLAT" & LOO_Weight > 0.5)) %>%
    ungroup() %>%
    filter(Model_Year %in% c("EXPAR1", "AR1")) %>%
    mutate(Viable = Delta_LOOIC < 10)

## ---------------- Summarising the models

# Computing the estimate summaries
mods <-
    mods %>%
    mutate(Summary = future_map(ID,
                                load_and_summarise,
                                progress = FALSE,
                                .progress = TRUE))

# Some formatting
est <-
    mods %>%
    select(DataID, Species, Population, Dist, Model_Year, Viable, Summary) %>%
    unnest(Summary) %>%
    mutate(Parameter = recode(Parameter,
                              t_intercept       = "theta",
                              log_w_intercept   = "log_Wmax"))

saveRDS(est, file = "../estimates.rds", version = 2)

## ---------------- Plotting the graphs

est <- readRDS("../estimates.rds")

## Format factors for labels
graph <-
    est %>%
    left_join(ids) %>%
    mutate(Parameter = if_else(Dist == "Binom",
                               recode(Parameter,
                                      log_Wmax = "logit_Wmax",
                                      log_w_sigma = "logit_w_sigma"),
                               Parameter),
           Parameter = factor(Parameter,
                              levels = c("theta", "t_sigma",
                                         "log_Wmax", "logit_Wmax",
                                         "log_w_sigma", "logit_w_sigma",
                                         "slope", "sl_sigma", "intercept", "int_sigma",
                                         "omega", "phi", "p_zi", "ind_sigma")),
           ID = fct_rev(ID))

levels(graph[["Parameter"]]) <-
    c("italic(θ)",
      "italic(σ[θ])",
      "log(italic(W)['max'])",
      "logit(italic(W)['max'])",
      "italic(σ)[log(italic(W)['max'])]",
      "italic(σ)[logit(italic(W)['max'])]",
      "italic(b)",
      "italic(σ[b])",
      "italic(a)",
      "italic(σ[a])",
      "italic(ω)",
      "italic(φ)",
      "italic(p)['zi']",
      "italic(σ)['ind']")
graph[["Dist"]] <-
    case_when(str_detect(graph[["Dist"]], "Poisson")    ~ "Poisson",
              str_detect(graph[["Dist"]], "Binom")      ~ "Binom")


## Creating the graphs
p_pois_dir <- 
    graph %>% filter(Dist == "Poisson", Model_Year == "EXPAR1") %>%
    {
        ggplot(.) +
            geom_linerange(data = ~ filter(.x, Viable),
                           aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_linerange(data =~  filter(.x, !Viable),
                           aes(x = ID, ymin = Low, ymax = Up), col = "grey70", size = 1) +
            geom_point(data = ~ filter(.x, Viable),
                       aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(data = ~ filter(.x, !Viable),
                       aes(x = ID, y = Mode), col = "#aa9696", size = 2) +
            geom_point(data = ~ filter(.x, Viable),
                       aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            geom_point(data = ~ filter(.x, !Viable),
                       aes(x = ID, y = Median), col = "#6d7494", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free_x",
                       nrow        = 3,
                       ncol        = 4,
                       labeller    = label_parsed) +
            coord_flip() +
            scale_x_discrete(limits = chuck(., "ID") %>%
                                      unique() %>%
                                      gtools::mixedsort()) + 
            ylab("Point estimates") + xlab("Dataset") +
            theme(text         = element_text(family = "Linux Biolinum O"),
                  axis.text.y  = element_text(size = 8),
                  axis.text.x  = element_text(size = 14),
                  axis.title   = element_text(size = 16),
                  strip.text   = element_text(family = "Linux Libertine O", size = 22))
    }

p_pois_opt <- 
    graph %>% filter(Dist == "Poisson", Model_Year == "AR1") %>%
    {
        ggplot(.) +
            geom_linerange(data = ~ filter(.x, Viable),
                           aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_linerange(data = ~ filter(.x, !Viable),
                           aes(x = ID, ymin = Low, ymax = Up), col = "grey70", size = 1) +
            geom_point(data = ~ filter(.x, Viable),
                       aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(data = ~ filter(.x, !Viable),
                       aes(x = ID, y = Mode), col = "#aa9696", size = 2) +
            geom_point(data = ~ filter(.x, Viable),
                       aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            geom_point(data = ~ filter(.x, !Viable),
                       aes(x = ID, y = Median), col = "#6d7494", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free_x",
                       nrow        = 3,
                       ncol        = 4,
                       labeller    = label_parsed) +
            coord_flip() +
            scale_x_discrete(limits = chuck(., "ID") %>%
                                      unique() %>%
                                      gtools::mixedsort()) + 
            ylab("Point estimates") + xlab("Dataset") +
            theme(text         = element_text(family = "Linux Biolinum O"),
                  axis.text.y  = element_text(size = 8),
                  axis.text.x  = element_text(size = 14),
                  axis.title   = element_text(size = 16),
                  strip.text   = element_text(family = "Linux Libertine O", size = 22))
    }

p_binom_dir <-
    graph %>% filter(Dist == "Binom", Model_Year == "EXPAR1") %>%
    {
        ggplot(.) +
            geom_linerange(data = ~ filter(.x, Viable),
                           aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_linerange(data =~  filter(.x, !Viable),
                           aes(x = ID, ymin = Low, ymax = Up), col = "grey70", size = 1) +
            geom_point(data = ~ filter(.x, Viable),
                       aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(data = ~ filter(.x, !Viable),
                       aes(x = ID, y = Mode), col = "#aa9696", size = 2) +
            geom_point(data = ~ filter(.x, Viable),
                       aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            geom_point(data = ~ filter(.x, !Viable),
                       aes(x = ID, y = Median), col = "#6d7494", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free_x",
                       ncol        = 4,
                       labeller    = label_parsed) +
            coord_flip() +
            scale_x_discrete(limits = chuck(., "ID") %>%
                                      unique() %>%
                                      gtools::mixedsort()) + 
            ylab("Point estimates") + xlab("Dataset") +
            theme(text         = element_text(family = "Linux Biolinum O"),
                  axis.text.y  = element_text(size = 8),
                  axis.text.x  = element_text(size = 14),
                  axis.title   = element_text(size = 16),
                  strip.text   = element_text(family = "Linux Libertine O", size = 22))
    }

p_binom_opt <-
    graph %>% filter(Dist == "Binom", Model_Year == "AR1") %>%
    {
        ggplot(.) +
            geom_linerange(data = ~ filter(.x, Viable),
                           aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_linerange(data =~  filter(.x, !Viable),
                           aes(x = ID, ymin = Low, ymax = Up), col = "grey70", size = 1) +
            geom_point(data = ~ filter(.x, Viable),
                       aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(data = ~ filter(.x, !Viable),
                       aes(x = ID, y = Mode), col = "#aa9696", size = 2) +
            geom_point(data = ~ filter(.x, Viable),
                       aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            geom_point(data = ~ filter(.x, !Viable),
                       aes(x = ID, y = Median), col = "#6d7494", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free_x",
                       ncol        = 4,
                       labeller    = label_parsed) +
            coord_flip() +
            scale_x_discrete(limits = chuck(., "ID") %>%
                                      unique() %>%
                                      gtools::mixedsort()) + 
            ylab("Point estimates") + xlab("Dataset") +
            theme(text         = element_text(family = "Linux Biolinum O"),
                  axis.text.y  = element_text(size = 8),
                  axis.text.x  = element_text(size = 14),
                  axis.title   = element_text(size = 16),
                  strip.text   = element_text(family = "Linux Libertine O", size = 22))
    }

## Saving the graphs
cairo_pdf("../../Figures/Estimates_Poisson.pdf", height = 10, width = 20)
(p_pois_opt + labs(tag = "FluctCorrOpt (Poisson)\n")) + (p_pois_dir + labs(tag = "FluctCorrDir (Poisson)\n")) &
    theme(plot.tag.position = "top", plot.tag = element_text(size = 20, face = "bold"))
dev.off()

cairo_pdf("../../Figures/Estimates_Binom.pdf", height = 8, width = 20)
(p_binom_opt + labs(tag = "FluctCorrOpt (Binomial)\n")) + (p_binom_dir + labs(tag = "FluctCorrDir (Binomial)\n")) &
    theme(plot.tag.position = "top", plot.tag = element_text(size = 20, face = "bold"))
dev.off()
