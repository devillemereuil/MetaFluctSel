library(tidyverse)
library(magrittr)
library(rstan)
library(progress)
library(furrr)
library(patchwork)
plan(multiprocess)
options(mc.cores = ifelse(availableCores() > 10,
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
taxa <- readRDS("../../Data/taxa.rds")

# List all mammals
ls_mammals <- 
    inner_join(ids, taxa %>% filter(Taxon == "Mammal")) %>%
    chuck("ID")


# Subsetting the two most complete models
mods <-
    mods %>%
    left_join(loo %>% select(DataID, Species, Population, Model_Year, Delta_LOOIC, LOO_Weight)) %>%
    group_by(DataID, Species, Population) %>%
    filter(!any(Model_Year == "FLAT" & LOO_Weight > 0.5)) %>%
    mutate(Weight_Opt = 
               LOO_Weight[Model_Year %in% c("IID", "AR1", "VARINT")] %>%
                   sum(),
           Weight_Dir = 
               LOO_Weight[Model_Year %in% c("EXPIID", "EXPAR1", "EXPVARINT")] %>%
                   sum(),
           Ratio_Opt_Dir = Weight_Opt / Weight_Dir,
           Keep_Theta = Ratio_Opt_Dir > 1,
           Viable = Delta_LOOIC < 10) %>%
    ungroup() %>%
    filter(Model_Year %in% c("EXPAR1", "AR1"))

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
    select(DataID, Species, Population, Dist, Model_Year, Viable, Keep_Theta, Summary) %>%
    unnest(Summary) %>%
    mutate(Parameter = recode(Parameter,
                              t_intercept       = "theta",
                              log_w_intercept   = "log_Wmax"))

## ------------------ Adding the unstandardised theta

sd_within <-
    readRDS("../STAN_mods.rds") %>%
    select(DataID, Species, Population, Data) %>%
    distinct() %>%
    mutate(SD_Within = map(Data, "SD_Within") %>%
                       map_dbl(unique)) %>%
    select(-Data)
    

unstd <-
    est %>%
    filter(Model_Year == "AR1", Parameter %in% c("theta", "t_sigma")) %>%
    left_join(sd_within) %>%
    mutate_at(vars(Median:Up), ~ . * SD_Within) %>%
    mutate(Parameter = recode(Parameter,
                              theta = "theta_unstd",
                              t_sigma = "t_sigma_unstd")) %>%
    select(-SD_Within)

est <- bind_rows(est, unstd)

saveRDS(est, file = "../estimates.rds", version = 2)

## ---------------- Plotting the graphs

est <- readRDS("../estimates.rds")

## Format factors for labels
graph <-
    est %>%
    left_join(ids) %>%
    filter(!str_detect(Parameter, "_unstd")) %>%
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

add_spacer <- function(limits) {
    cut <- max(which(limits %in% ls_mammals))
    c(limits[1:cut], "", limits[(cut+1):length(limits)])
}


## Creating the graphs
p_pois_dir <- 
    graph %>% 
    filter(Dist == "Poisson", Model_Year == "EXPAR1") %>%
    mutate(ID = fct_drop(ID)) %>%
    {
        ggplot(.) +
            geom_linerange(data = ~ filter(.x, !Keep_Theta),
                           aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_linerange(data =~  filter(.x, Keep_Theta),
                           aes(x = ID, ymin = Low, ymax = Up), col = "grey70", size = 1) +
            geom_point(data = ~ filter(.x, !Keep_Theta),
                       aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(data = ~ filter(.x, Keep_Theta),
                       aes(x = ID, y = Mode), col = "#aa9696", size = 2) +
            geom_point(data = ~ filter(.x, !Keep_Theta),
                       aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            geom_point(data = ~ filter(.x, Keep_Theta),
                       aes(x = ID, y = Median), col = "#6d7494", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free_x",
                       nrow        = 3,
                       ncol        = 4,
                       labeller    = label_parsed) +
            coord_flip() +
            scale_x_discrete(breaks = levels(.[["ID"]]),
                             limits = add_spacer(levels(.[["ID"]]))) + 
            ylab("Point estimates") + xlab("Dataset") +
            theme(text         = element_text(family = "Linux Biolinum O"),
                  axis.text.y  = element_text(size = 8),
                  axis.text.x  = element_text(size = 14),
                  axis.title   = element_text(size = 16),
                  strip.text   = element_text(family = "Linux Libertine O", size = 22))
    }

p_pois_opt <- 
    graph %>% 
    filter(Dist == "Poisson", Model_Year == "AR1") %>%
    mutate(ID = fct_drop(ID)) %>%
    {
        ggplot(.) +
            geom_linerange(data = ~ filter(.x, Keep_Theta),
                           aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_linerange(data = ~ filter(.x, !Keep_Theta),
                           aes(x = ID, ymin = Low, ymax = Up), col = "grey70", size = 1) +
            geom_point(data = ~ filter(.x, Keep_Theta),
                       aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(data = ~ filter(.x, !Keep_Theta),
                       aes(x = ID, y = Mode), col = "#aa9696", size = 2) +
            geom_point(data = ~ filter(.x, Keep_Theta),
                       aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            geom_point(data = ~ filter(.x, !Keep_Theta),
                       aes(x = ID, y = Median), col = "#6d7494", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free_x",
                       nrow        = 3,
                       ncol        = 4,
                       labeller    = label_parsed) +
            coord_flip() +
            scale_x_discrete(breaks = levels(.[["ID"]]),
                             limits = add_spacer(levels(.[["ID"]]))) + 
            ylab("Point estimates") + xlab("Dataset") +
            theme(text         = element_text(family = "Linux Biolinum O"),
                  axis.text.y  = element_text(size = 8),
                  axis.text.x  = element_text(size = 14),
                  axis.title   = element_text(size = 16),
                  strip.text   = element_text(family = "Linux Libertine O", size = 22))
    }

p_binom_dir <-
    graph %>%
    filter(Dist == "Binom", Model_Year == "EXPAR1") %>%
    mutate(ID = fct_drop(ID)) %>%
    {
        ggplot(.) +
            geom_linerange(data = ~ filter(.x, !Keep_Theta),
                           aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_linerange(data =~  filter(.x, Keep_Theta),
                           aes(x = ID, ymin = Low, ymax = Up), col = "grey70", size = 1) +
            geom_point(data = ~ filter(.x, !Keep_Theta),
                       aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(data = ~ filter(.x, Keep_Theta),
                       aes(x = ID, y = Mode), col = "#aa9696", size = 2) +
            geom_point(data = ~ filter(.x, !Keep_Theta),
                       aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            geom_point(data = ~ filter(.x, Keep_Theta),
                       aes(x = ID, y = Median), col = "#6d7494", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free_x",
                       ncol        = 4,
                       labeller    = label_parsed) +
            coord_flip() +
            scale_x_discrete(breaks = levels(.[["ID"]]),
                             limits = levels(.[["ID"]])) + 
            ylab("Point estimates") + xlab("Dataset") +
            theme(text         = element_text(family = "Linux Biolinum O"),
                  axis.text.y  = element_text(size = 8),
                  axis.text.x  = element_text(size = 14),
                  axis.title   = element_text(size = 16),
                  strip.text   = element_text(family = "Linux Libertine O", size = 22))
    }

p_binom_opt <-
    graph %>% 
    filter(Dist == "Binom", Model_Year == "AR1") %>%
    mutate(ID = fct_drop(ID)) %>%
    {
        ggplot(.) +
            geom_linerange(data = ~ filter(.x, Keep_Theta),
                           aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_linerange(data =~  filter(.x, !Keep_Theta),
                           aes(x = ID, ymin = Low, ymax = Up), col = "grey70", size = 1) +
            geom_point(data = ~ filter(.x, Keep_Theta),
                       aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(data = ~ filter(.x, !Keep_Theta),
                       aes(x = ID, y = Mode), col = "#aa9696", size = 2) +
            geom_point(data = ~ filter(.x, Keep_Theta),
                       aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            geom_point(data = ~ filter(.x, !Keep_Theta),
                       aes(x = ID, y = Median), col = "#6d7494", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free_x",
                       ncol        = 4,
                       labeller    = label_parsed) +
            coord_flip() +
            scale_x_discrete(breaks = levels(.[["ID"]]),
                             limits = levels(.[["ID"]])) + 
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

## ------------------------- Plotting the unstandardised estimates

graph_unstd <-
    est %>%
    filter(str_detect(Parameter, "unstd")) %>%
    left_join(ids) %>%
    group_by(ID) %>%
    mutate(LargeSD = Median[Parameter == "t_sigma_unstd"] > 50) %>%
    ungroup() %>%
    mutate(Parameter = factor(Parameter,
                              levels = c("theta_unstd", "t_sigma_unstd")) %>%
                       recode(theta_unstd = "italic(θ[paste('unstd')])",
                              t_sigma_unstd = "italic(σ[θ[paste('unstd')]])"),
           ID = fct_rev(ID)) %>%
    filter(Keep_Theta)

make_graph <-
    . %>% 
    { 
        ggplot(.) +
            geom_linerange(aes(x = ID, ymin = Low, ymax = Up), col = "black", size = 1) +
            geom_point(aes(x = ID, y = Mode), col = "#aa0000", size = 2) +
            geom_point(aes(x = ID, y = Median), col = "#0055ff", size = 2) +
            facet_wrap(~ Parameter,
                       scale       = "free",
                       ncol        = 2,
                       labeller    = label_parsed) +
            coord_flip() +
            ylab("Point estimates") + xlab("Dataset") +
            theme(text         = element_text(family = "Linux Biolinum O"),
                  axis.text.y  = element_text(size = 8),
                  axis.text.x  = element_text(size = 14),
                  axis.title   = element_text(size = 16),
                  strip.text   = element_text(family = "Linux Libertine O", size = 22),
                  strip.text.y = element_blank(),
                  strip.background.y = element_blank())
    }

cairo_pdf("../../Figures/Unstandardised_Theta.pdf", height = 6, width = 8)
    (
        make_graph(graph_unstd %>% filter(!LargeSD)) +
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_text(hjust = 0.4))
    ) / (
        make_graph(graph_unstd %>% filter(LargeSD)) +
            xlab("") +
            theme(strip.background = element_blank(),
                  strip.text = element_blank())
    ) + plot_layout(heights = c(sum(!graph_unstd[["LargeSD"]]), sum(graph_unstd[["LargeSD"]])))
dev.off()
