library(tidyverse)
library(magrittr)
library(rstan)
library(progress)
library(furrr)
library(brms)
plan(multiprocess)
options(mc.cores = ifelse(availableCores() > 10,
                          10,
                          availableCores()))
filter  <- dplyr::filter
lag     <- dplyr::lag
extract <- rstan::extract

theme_set(theme_bw())

source("00_functions.R")

## ---------------------------------- Loading the data

# Loading the meta-data
mods <- readRDS("../STAN_mods.rds")
# Loading the LOO scores
loo  <- readRDS("../loo.rds")
# Getting the IDs
ids <- readRDS("../FinalIDs.rds")
# Getting the taxon information
taxa <- readRDS("../../Data/taxa.rds")
# Getting the gradients
gradients <- readRDS("../gradients.rds")

# List all mammal species (needed for spacer)
list_birds <- 
    ids %>% 
    left_join(taxa) %>%
    filter(Taxon == "Bird") %>%
    chuck("ID") %>%
    as.character()
list_mammals <- 
    ids %>% 
    left_join(taxa) %>% 
    filter(Taxon == "Mammal") %>% 
    chuck("ID") %>%
    as.character()

# Selecting viable optimum models
mods <-
    loo %>% 
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
    select(-LOOIC, -Delta_LOOIC, -Weight_Opt, -Weight_Fluct, -Ratio_Opt_Fluct) %>%
    inner_join(mods)

# Getting the thetas (opt models) and beta (dir models) with fluctuations
track <-
    mods %>%
    filter(Model_Year %in% "AR1") %>%
    select(DataID, Species, Population, Data, Viable, Keep_Theta) %>%
    mutate(Theta = future_map(mods[["ID"]][mods[["Model_Year"]] == "AR1"],
                              load_and_get_theta,
                              progress = FALSE,
                              .progress = TRUE),
           Beta  = future_map(mods[["ID"]][mods[["Model_Year"]] == "EXPAR1"],
                              load_and_get_estexp,
                              progress = FALSE,
                              .progress = TRUE) %>%
                   map("Slope")) %>%
    filter(Viable) %>%
    select(-Viable)

## ---------------------------------- Computing the deltas

# Using the phenotypic data
track <-
    track %>%
    mutate(Deltas = map(Data, compute_delta_zbar, column = "Pheno_scale"))

# Using thetas
track <-
    track %>%
    mutate(Delta_Theta = map(Theta, ~ {.[ , 2:ncol(.)] - .[ , 1:(ncol(.) - 1)]}),
           Delta_Beta  = map(Beta,  ~ {.[ , 2:ncol(.)] - .[ , 1:(ncol(.) - 1)]}))

# Cleaning up a few columns and removing datasets with too little information
track <-
    track %>%
    left_join(ids) %>%
    select(-DataID, -Species, -Population, -Data, -Theta, -Beta) %>%
    mutate(ToKeep = map(Deltas, ~ which(.[["N_CommonID"]] >= 5)),
           Deltas = map2(Deltas, ToKeep, ~ slice(.x, .y)),
           Delta_Theta = map2(Delta_Theta, ToKeep, ~ .x[ , .y]),
           Delta_Beta = map2(Delta_Beta, ToKeep, ~ .x[ , .y]),
           ToKeep = NULL) %>%
    filter(map_int(Deltas, nrow) >= 10)

# Checking that the subsetting is not too off!
track %$% 
    map_dbl(Deltas, ~ cor(.[["Zbar"]], .[["CommonZbar"]])) %T>%
    {print(any(. < 0.80))}

## ---------------------------------- Comparing the deltas

# Compute the correlations between delta_theta and delta_zbar
nitt <- 4000
track <-
    track %>%
    # We need to Monte Carlo sample the Zbar to account for uncertainty in mean
    # NOTE This assumes a Gaussian distribution of the Zbar...
    mutate(DeltaZbar_Dist = future_map(Deltas,
                                       # Rerun will produce nitt items in a list
                                       ~ rerun(nitt,
                                               rnorm(nrow(.),
                                                     mean = .[["DeltaCommonZbar"]],
                                                     sd = .[["DeltaCommonZbar_SE"]])) %>%
                                         # Now we bind the nitt vectors into a matrix
                                         {do.call("rbind", .)}),
           # Now we can compute the correlation row by row between the two matrices
           Corr_Theta = future_map2(DeltaZbar_Dist, Delta_Theta,
                                rowwise_correlate,
                                method = "pearson"),
           Corr_Beta = future_map2(DeltaZbar_Dist, Delta_Beta,
                                rowwise_correlate,
                                method = "pearson"),
           Pval_Theta = future_map_dbl(Corr_Theta, compute_pval),
           Pval_Beta = future_map_dbl(Corr_Beta, compute_pval))

saveRDS(track, file = "../tracking.rds", version = 2)

## ---------------------------------- Now plotting the results

add_spacer <- function(limits) {
    cut <- max(which(limits%in% list_mammals))
    c(limits[1:cut], "", limits[(cut+1):length(limits)])
}

## Plotting the correlations
# Setting up the data
graph_corr <-
    track %>%
    filter(Keep_Theta) %>%
    mutate(ID = factor(ID, levels =  c(list_birds, list_mammals)) %>%
                fct_rev() %>%
                fct_drop()) %>%
    select(ID, Corr_Theta, Pval_Theta, Corr_Beta, Pval_Beta) %>%
    unnest_legacy()

# Creating the plot
p_corr_theta <-
    ggplot(graph_corr) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_violin(aes(x = ID, y = Corr_Theta, fill = Pval_Theta < 0.05)) +
    coord_flip() +
    annotate("text",
             y = -1,
             x = n_distinct(graph_corr[["ID"]]) + 1.5,
             label = "Birds",
             face = "bold",
             size = 8,
             hjust= 0,
             family = "Linux Biolinum O") +
    annotate("text",
             y = -1,
             x = max(which(levels(graph_corr[["ID"]]) %in% list_mammals)) + 0.75,
             label = "Mammals",
             face = "bold",
             size = 8,
             hjust = 0,
             family = "Linux Biolinum O") +
    ylab("Plasticity - optimum correlation") + xlab("Dataset") +
    scale_x_discrete(breaks = levels(graph_corr[["ID"]]),
                     limits = add_spacer,
                     expand = expansion(add = 1)) +
    scale_fill_manual(values = c("#00557f", "#C65A68"), name = "Corr. ≠ 0?") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 10),
          axis.text.x  = element_text(size = 18),
          axis.title   = element_text(size = 20),
          panel.grid.major.y = element_blank(),
          strip.text   = element_text(family = "Linux Biolinum O", size = 22))


# Now saving the graphs
cairo_pdf("../../Figures/Tracking_correlations_theta.pdf", height = 8, width = 6)
plot(p_corr_theta + theme(legend.position = "none"))
dev.off()

# Creating the plot
p_corr_beta <-
    ggplot(graph_corr) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_violin(aes(x = ID, y = Corr_Beta, fill = Pval_Beta < 0.05)) +
    coord_flip() +
    ylab("Plasticity - gradient correlation") + xlab("Dataset") +
    scale_fill_manual(values = c("#00557f", "#C65A68"), name = "Corr. ≠ 0?") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 10),
          axis.text.x  = element_text(size = 18),
          axis.title   = element_text(size = 20),
          strip.text   = element_text(family = "Linux Biolinum O", size = 22))


# Now saving the graphs
cairo_pdf("../../Figures/Tracking_correlations_beta.pdf", height = 8, width = 6)
plot(p_corr_beta + theme(legend.position = "none"))
dev.off()

p_corr_talk <-
    ggplot(graph_corr) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_violin(aes(x = ID, y = Corr_Theta, fill = Pval_Theta < 0.05)) +
    coord_flip() +
    annotate("text",
             y = -1,
             x = n_distinct(graph_corr[["ID"]]) + 1.5,
             label = "Birds",
             fontface = "bold",
             size = 8,
             hjust= 0,
             family = "Noto Sans") +
    annotate("text",
             y = -1,
             x = max(which(levels(graph_corr[["ID"]]) %in% list_mammals)) + 0.75,
             label = "Mammals",
             fontface = "bold",
             size = 8,
             hjust = 0,
             family = "Noto Sans") +
    ylab("Plasticity - optimum correlation") + xlab("Dataset") +
    scale_fill_manual(values = c("#00557f", "#AA0000"), name = "Corr. ≠ 0?") +
    scale_x_discrete(breaks = levels(graph_corr[["ID"]]),
                     limits = add_spacer,
                     expand = expansion(add = 1)) +
    theme(text = element_text(family = "Noto Sans", size = 16),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())
    

cairo_pdf("../../Figures/Tracking_correlations_preso.pdf", height = 8, width = 5)
plot(p_corr_talk + theme(legend.position = "none"))
dev.off()

## Plotting the actual relationship

# Setting up the data
graph_est <-
    track %>%
    mutate(Delta_Zbar = map(Deltas, "DeltaCommonZbar"),
           Delta_Theta = map(Delta_Theta, ~ apply(., 2, median)),
           ID = as.factor(ID) %>% fct_rev()) %>%
    select(ID, Delta_Zbar, Delta_Theta) %>%
    unnest_legacy()

# Creating the plot
p_est <-
    ggplot(graph_est) +
    geom_point(aes(x = Delta_Theta, y = Delta_Zbar)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey", alpha = 0.5) +
    facet_wrap(~ ID, scale = "free") +
    xlab(bquote(Delta~theta)) + ylab(bquote(Delta~bar(z))) +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 8),
          axis.text.x  = element_text(size = 8),
          axis.title   = element_text(family = "Linux Libertine O", size = 16))

# Now saving the graphs
cairo_pdf("../../Figures/Tracking_estimates.pdf", height = 12, width = 10)
plot(p_est)
dev.off()

## ---------------------------------- Some statistical tests

## Formatting the MCMC for "multiple imputation"
df_corr <-
    track %>%
    filter(Keep_Theta) %>%
    mutate(ID = as.factor(ID) %>% fct_rev(),
           Sample_Size = map_int(Deltas, nrow)) %>%
    left_join(ids %>% select(-DataID)) %>%
    left_join(taxa) %>%
    select(ID, Species, Population, Taxon, Corr_Theta, Sample_Size, Pval_Theta) %>%
    unnest_legacy() %>%
    group_by(ID) %>%
    sample_n(100 * Sample_Size, replace = TRUE) %>%
    mutate(Sample = rep(1:100, each = unique(Sample_Size))) %>%
    ungroup()

dfs_tot <-
    df_corr %>%
    group_by(Sample) %>%
    nest_legacy(.key = "Df") %>%
    pluck("Df")

## Running multiple models for all datasets
prior <- c(
    prior(normal(0, 1), "b"),
    prior(normal(0, 1), "sd"),
    prior(normal(0, 1), "sigma")
)
form <- brmsformula(Corr_Theta ~ 1 + Taxon + (1|ID),
                    sparse = TRUE)
model_tot <- 
    brm_multiple(formula    = form,
                 data       = dfs_tot,
                 save_ranef = FALSE,
                 chains     = 2,
                 iter       = 3000,
                 warmup     = 2000,
                 thin       = 50,
                 prior      = prior)
saveRDS(model_tot, file = "../model_corr_tot.rds", version = 2)


## Compute the p-values
model_tot <- readRDS("../model_corr_tot.rds")

# For all datasets
as.data.frame(model_tot) %>%
    mutate(mu_Bird      = b_Intercept,
           mu_Mammal    = b_Intercept + b_TaxonMammal) %>%
    summarise_at(vars(contains("mu")),
                 list(median = median,
                      CI = ~ coda::HPDinterval(as.mcmc(.)),
                      pval = compute_pval)) %>%
    pivot_longer(everything(),
                 names_pattern = "mu_(.+)_(.+)",
                 names_to = c("Taxon", ".value"),
                 values_to = "Value")
# Taxon  median  CI[,1]  [,2]   pval
# Bird    0.247  0.0722 0.437 0.0095
# Mammal  0.133 -0.168  0.426 0.350 

