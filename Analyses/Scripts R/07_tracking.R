library(tidyverse)
library(magrittr)
library(rstan)
library(progress)
library(furrr)
library(brms)
plan(multiprocess)
options(mc.cores = ifelse(availableCores() > 4,
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
# Selecting viable optimum models
mods <- inner_join(mods,
                   loo %>% 
                       group_by(DataID, Species, Population) %>%
                       filter(!any(Model_Year == "FLAT" & LOO_Weight > 0.5)) %>%
                       ungroup() %>%
                       filter(Delta_LOOIC < 10) %>% 
                       select(-LOOIC, -Delta_LOOIC))

# Getting the thetas for the opt models with fluctuations
track <-
    mods %>%
    filter(Model_Year %in% c("AR1")) %>%
    mutate(Theta = future_map(ID,
                              load_and_get_theta,
                              progress = FALSE,
                              .progress = TRUE))

## ---------------------------------- Computing the deltas

# Using the phenotypic data
track <-
    track %>%
    mutate(Deltas = map(Data, compute_delta_zbar, column = "Pheno_scale"))

# Using thetas
track <-
    track %>%
    mutate(Delta_Theta = map(Theta, ~ {.[ , 2:ncol(.)] - .[ , 1:(ncol(.) - 1)]}))

# Cleaning up a few columns and removing datasets with too little information
track <-
    track %>%
    select(-ID) %>%
    left_join(ids) %>%
    select(-DataID, -Species, -Population, -Data, -Model_Path, -Theta) %>%
    mutate(ToKeep = map(Deltas, ~ which(.[["N_CommonID"]] >= 5)),
           Deltas = map2(Deltas, ToKeep, ~ slice(.x, .y)),
           Delta_Theta = map2(Delta_Theta, ToKeep, ~ .x[ , .y]),
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
           Corr_P = future_map2(DeltaZbar_Dist, Delta_Theta,
                                rowwise_correlate,
                                method = "pearson"),
           Corr_S = future_map2(DeltaZbar_Dist, Delta_Theta,
                                rowwise_correlate,
                                method = "spearman"),
           Pval = future_map_dbl(Corr_P, compute_pval))

saveRDS(track, file = "../tracking.rds", version = 2)

## ---------------------------------- Now plotting the results

## Plotting the correlations
# Setting up the data
graph_corr <-
    track %>% 
    mutate(ID = as.factor(ID) %>% fct_rev()) %>%
    select(ID, Corr_P, Pval) %>%
    unnest_legacy()

# Creating the plot
p_corr <-
    ggplot(graph_corr) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_violin(aes(x = ID, y = Corr_P, fill = Pval < 0.05)) +
    coord_flip() +
    ylab("Plasticity - optimum correlation") + xlab("Dataset") +
    scale_fill_manual(values = c("#00557f", "#C65A68"), name = "Corr. ≠ 0?") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 10),
          axis.text.x  = element_text(size = 18),
          axis.title   = element_text(size = 20),
          strip.text   = element_text(family = "Linux Biolinum O", size = 22))


# Now saving the graphs
cairo_pdf("../../Figures/Tracking_correlations.pdf", height = 8, width = 6)
plot(p_corr + theme(legend.position = "none"))
dev.off()

p_corr_talk <-
    ggplot(graph_corr) +
    geom_violin(aes(x = ID, y = Corr_P, fill = Pval < 0.05)) +
    coord_flip() +
    ylab("Plasticity - optimum correlation") + xlab("Dataset") +
    scale_fill_manual(values = c("#00557f", "#AA0000"), name = "Corr. ≠ 0?") +
    theme(text = element_text(family = "Noto Sans", size = 16),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
# ,
#           axis.text.y  = element_text(size = 8),
#           axis.text.x  = element_text(size = 14),
#           axis.title   = element_text(size = 16))

cairo_pdf("../../Figures/Tracking_correlations_preso.pdf", height = 10, width = 5)
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
    mutate(ID = as.factor(ID) %>% fct_rev(),
           Sample_Size = map_int(Deltas, nrow)) %>%
    left_join(ids %>% select(-DataID)) %>%
    left_join(taxa) %>%
    select(ID, Species, Population, Taxon, Corr_P, Sample_Size, Pval) %>%
    unnest_legacy() %>%
    group_by(ID) %>%
    sample_n(100 * Sample_Size, replace = TRUE) %>%
    mutate(Sample = rep(1:100, each = unique(Sample_Size))) %>%
    ungroup()

# Using all the datasets
dfs_tot <-
    df_corr %>%
    group_by(Sample) %>%
    nest_legacy(.key = "Df") %>%
    pluck("Df")

# Using the non-significant datasets
dfs_sub <-
    df_corr %>%
    filter(Pval > 0.05) %>%
    group_by(Sample) %>%
    nest_legacy(.key = "Df") %>%
    pluck("Df")

## Running multiple models for all datasets
prior <- c(
    prior(normal(0, 1), "b"),
    prior(normal(0, 1), "sd"),
    prior(normal(0, 1), "sigma")
)
form <- brmsformula(Corr_P ~ 1 + Taxon + (1|ID),
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

# Running multiple models for non-significant datasets
model_sub <- update(model_tot, newdata = dfs_sub)
saveRDS(model_sub, file = "../model_corr_sub.rds", version = 2)

## Compute the p-values
model_tot <- readRDS("../model_corr_tot.rds")
model_sub <- readRDS("../model_corr_sub.rds")
# For all datasets
as.data.frame(model_tot) %>%
    mutate(mu_Bird      = b_Intercept,
           mu_Mammal    = b_Intercept + b_TaxonMammal) %>%
    summarise_at(vars(contains("mu")),
                 list(median = median,
                      pval = compute_pval))
#  mu_Bird_median mu_Mammal_median mu_Bird_pval mu_Mammal_pval
#       0.2077379        0.1146359        5e-04         0.2755

# For non-significant datasets
as.data.frame(model_sub) %>%
    mutate(mu_Bird      = b_Intercept,
           mu_Mammal    = b_Intercept + b_TaxonMammal) %>%
    summarise_at(vars(contains("mu")),
                 list(median = median,
                      pval = compute_pval))
#  mu_Bird_median mu_Mammal_median mu_Bird_pval mu_Mammal_pval
#1      0.1140279        0.1141573       0.0215          0.166
