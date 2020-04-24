library(tidyverse)
library(magrittr)
library(rstan)
library(brms)
library(progress)
library(furrr)
library(ggrepel)
library(ggnewscale)
library(grImport)
library(patchwork)
plan(multiprocess)
options(mc.cores = ifelse(availableCores() > 4,
                          10,
                          availableCores()))
filter  <- dplyr::filter
lag     <- dplyr::lag
extract <- rstan::extract

theme_set(theme_bw())
rstan_options(auto_write = TRUE)

source("00_functions.R")

## ---------------------------------- Loading the data
# Loading the meta-data
mods <- readRDS("../STAN_mods.rds")
# Loading the LOO scores
loo  <- readRDS("../loo.rds")
# Loading the IDs
ids  <- readRDS("../FinalIDs.rds")
# Loading the taxa
taxa <- readRDS("../../Data/taxa.rds")


# Subsetting the two most complete models
mods <-
    mods %>%
    left_join(loo %>% select(DataID, Species, Population, Model_Year, Delta_LOOIC, LOO_Weight)) %>%
    group_by(DataID, Species, Population) %>%
    filter(!any(Model_Year == "FLAT" & LOO_Weight > 0.5)) %>%
    ungroup() %>%
    filter(Model_Year %in% c("EXPAR1", "AR1")) %>%
    mutate(Viable = Delta_LOOIC < 10)

print("Handling optimum models")
# Obtaining the optimum locations from the models for which there is an optimum
opt <-
    mods %>%
    filter(Model_Year == "AR1") %>%
    mutate(Theta = future_map(ID,
                              load_and_get_theta,
                              progress = FALSE,
                              .progress = TRUE),
           Theta = map_if(Theta,
                          ~ length(dim(.)) == 1,
                          ~ matrix(., ncol = 1)),
           Omega = future_map(ID,
                              load_and_get_omega,
                              progress = FALSE,
                              .progress = TRUE),
           SD_Theta = future_map(ID,
                                 load_and_get_t_sigma,
                                 progress = FALSE,
                                 .progress = TRUE))

print("Handling directional models")
# Obtaining the gradients (slope) for the directional models
dir <-
    mods %>%
    filter(Model_Year == "EXPAR1") %>%
    mutate(EstEXP = future_map(ID,
                               load_and_get_estexp,
                               progress = FALSE,
                               .progress = TRUE)) %>% 
    unnest_wider(EstEXP) %>%
    mutate(Nyear    = map_int(Data, ~ length(full_seq(.[["Year"]], 1))),  
           Intercept = map_if(Intercept,
                              ~ length(dim(.)) == 1,
                              ~ matrix(., ncol = 1)) %>%
                       # Ensure formatting if ever using EXPNULL
                       map2(Nyear,
                            ~ {matrix(.x, ncol = .y, nrow = nrow(.x))}),
           Slope     = map_if(Slope,
                              ~ length(dim(.)) == 1,
                              ~ matrix(., ncol = 1)) %>%
                       # Ensure formatting if ever using EXPNULL
                       map2(Nyear,
                            ~ {matrix(.x, ncol = .y, nrow = nrow(.x))}),
           Nyear    = NULL)

## ---------------------------------- Computing the mismatches

print("Computing the mistmaches")
# Computing the average lay date (zbar), mismatch and gradient (beta) for each year
opt <-
    opt %>%
    mutate(Zbar     = map(Data, ~ compute_zbar(., "Pheno_scale")),
           SD_Zbar  = map_dbl(Zbar, ~ sd(., na.rm = TRUE)),
           CorrTZ   = future_map2(Theta, Zbar,
                                  ~ apply(array(c(.x, .y), dim = c(dim(.x), 2)),
                                          1,
                                          function(mat) {cor(mat[ , 1], mat[ , 2], use = "complete.obs")})),
           Mismatch = map2(Theta, Zbar,
                           ~ {matrix(.x, ncol = length(.y), nrow = nrow(.x)) -
                              matrix(.y, ncol = length(.y), nrow = nrow(.x), byrow = TRUE)}),
           Beta     = map2(Mismatch, Omega,
                           ~ {.x / (matrix(.y, ncol = ncol(.x), nrow = nrow(.x))^2 + 1)}),
           Beta_Nopl= map2(Theta, Omega,
                           ~ {.x / (matrix(.y, ncol = ncol(.x), nrow = nrow(.x))^2 + 1)}),
           Var_Beta = map(Beta, ~ apply(., 1, var, na.rm = TRUE)),
           Var_Beta_Nopl    = map(Beta_Nopl, ~ apply(., 1, var, na.rm = TRUE)),
           SD_Beta          = map(Beta, ~ apply(., 1, sd, na.rm = TRUE)),
           SD_Beta_Nofl     = map2(SD_Zbar, Omega,
                                   ~ as.vector((.x) / ((.y)^2 + 1))),
           SD_Beta_Nopl     = map(Beta_Nopl, ~ apply(., 1, sd, na.rm = TRUE)),
           SD_Beta_Notr     = pmap(list(SD_Theta, SD_Zbar, Omega),
                                   ~ as.vector(sqrt((..1)^2 + (..2)^2) / ((..3)^2 + 1))),
           SD_Beta_Hypo     = pmap(list(SD_Theta, SD_Zbar, CorrTZ, Omega),
                                   ~ as.vector(sqrt((..1)^2 + (..2)^2 - (2 * ..3 * ..1 * ..2)) / ((..4)^2 + 1))))

# Formatting gradients for directional models
dir <-
    dir %>%
    # Complex way of compute beta because we need to account for 
    mutate(Pheno  = map(Data,
                        ~ select(., Year, Pheno_scale) %>% 
                          nest(z = Pheno_scale) %>%
                          complete(Year = full_seq(Year, 1)) %>%
                          transpose() %>%
                          map(c("z", "Pheno_scale"))),
           Estimates_list = map2(Intercept, Slope,
                                 ~ array(c(.x, .y), dim = c(dim(.x), 2)) %>%
                                   apply(2, list) %>%
                                   flatten()),
           Beta     = if_else(Dist == "Binom",
                              # Computing for a logit, needs to use the phenotypic data
                              # Beta = Slope * (1 - mean(W(z)²) / mean(W(z))) see compute_beta_logit
                              future_map2(Estimates_list, Pheno,
                                          ~ map2_dfc(.x, .y,
                                                 function(arg1, arg2) {
                                                     apply(arg1, 1, compute_beta_logit, z = arg2)
                                                 }) %>%
                                            as.matrix(),
                                          .progress = TRUE),
                              # For other Poisson-based models, the Slope is already Beta...
                              future_map(Slope, ~ .x)),
           Estimates_list = NULL,
           Pheno  = NULL,
           # Computing the rest is easy
           Var_Beta = map(Beta, ~ apply(., 1, var, na.rm = TRUE)),
           SD_Beta  = map(Beta, ~ apply(., 1, sd, na.rm = TRUE)))

# Merging all gradients
all_grad <-
    bind_rows(opt %>% select(-Theta, -Omega, -Mismatch),
              dir %>% select(-Intercept, -Slope)) %>%
    mutate(Var_Sign = future_map(Beta, ~ apply(., 1, compute_var_sign)))

# betas <- map(opt[["Beta"]], ~ apply(., 2, median))
# thetas <- map(opt[["Theta"]], ~ apply(., 2, median))
# zbars <- opt[["Zbar"]]
# test <- 
#     mods %>% 
#     filter(Model_Year == "AR1") %>% 
#     mutate(CorrBT = map2_dbl(betas, thetas, ~ cor(.x, .y, use = "complete.obs")),
#            CorrBZ = map2_dbl(betas, zbars, ~ cor(.x, .y, use = "complete.obs")))
# 
# test %>% 
#     filter(CorrBT < 0.8) %>% 
#     select(DataID, Population, Species, Viable, CorrBT) %>%
#     print(n = Inf)
# 
# test %>% 
#     filter(CorrBZ > -0.8) %>% 
#     select(DataID, Population, Species, Viable, CorrBZ) %>%
#     print(n = Inf)
# 
# test %>% 
#     filter(CorrBT < 0.8, CorrBZ > -0.8) %>% 
#     select(DataID, Population, Species, Viable, CorrBT, CorrBZ)

## ---------------------------------- Formatting the posterior distributions

print("Formatting the posterior distributions")
# Function to compute point estimates and format them from a matrix
compute_pe_matrix <-
    . %>%
    apply(2, compute_point_estimates) %>%
    transpose() %>%
    map(unlist) %>%
    {tibble(!!!.)}

# Extracting point estimates
point_grad <-
    all_grad %>%
    mutate(Mean_Beta    = map(Beta, ~ apply(., 1, mean, na.rm = TRUE)),
           Beta         = map(Beta, compute_pe_matrix),
           CV_Beta      = map2(Mean_Beta, Var_Beta, ~ (sqrt(.y)/(.x))*100),
           Mean_Beta    = map(Mean_Beta, compute_point_estimates),
           Var_Beta     = map(Var_Beta, compute_point_estimates),
           SD_Beta      = map(SD_Beta, compute_point_estimates),
           SD_Beta_Nofl = map_if(SD_Beta_Nofl,
                                 ~ !is.null(.),
                                 compute_point_estimates,
                                 .else = ~ list(Median = NA,
                                                Mode = NA,
                                                SE = NA,
                                                Low = NA,
                                                Up = NA)),
           SD_Beta_Nopl = map_if(SD_Beta_Nopl,
                                 ~ !is.null(.),
                                 compute_point_estimates,
                                 .else = ~ list(Median = NA,
                                                Mode = NA,
                                                SE = NA,
                                                Low = NA,
                                                Up = NA)),
           SD_Beta_Notr = map_if(SD_Beta_Notr,
                                 ~ !is.null(.),
                                 compute_point_estimates,
                                 .else = ~ list(Median = NA,
                                                Mode = NA,
                                                SE = NA,
                                                Low = NA,
                                                Up = NA)),
           SD_Beta_Hypo = map_if(SD_Beta_Hypo,
                                 ~ !is.null(.),
                                 compute_point_estimates,
                                 .else = ~ list(Median = NA,
                                                Mode = NA,
                                                SE = NA,
                                                Low = NA,
                                                Up = NA)),
           CV_Beta      = map(CV_Beta, compute_point_estimates),
           Var_Sign     = map(Var_Sign, compute_point_estimates))

saveRDS(point_grad, file = "../gradients.rds", version = 2)

## ---------------------------------- Performing a sort of meta-analysis

print("Running the models")

## Formatting the MCMC for "multiple imputation"
# Getting the relevant MCMC samples for the "meta-analysis"
mcmc_meta <-
    opt %>%
    select(-ID) %>%
    left_join(ids) %>%
    left_join(taxa) %>%
    mutate(Mean_Beta    = map(Beta,  ~ apply(., 1, mean, na.rm = TRUE)),
           Mean_Theta   = map(Theta, ~ apply(., 1, mean, na.rm = TRUE)),
           Dist         = ifelse(str_detect(Dist, "Binom"), "Binom", "Poisson") %>%
                          factor(levels = c("Poisson", "Binom"))) %>%
    select(ID, Species, Taxon, Population, Dist, Mean_Theta, SD_Theta, Mean_Beta, SD_Beta, Omega) 

# Generating 100 datasets based on the MCMC samples
dfs_meta <-
    mcmc_meta %>%
    unnest_legacy() %>%
    group_by(ID) %>%
    sample_n(100) %>%
    mutate(Sample = 1:n()) %>%
    ungroup() %>%
    group_by(Sample) %>%
    nest_legacy(.key = "Df") %>%
    pluck("Df")

## Now running the five models
prior_n <- c(
    prior(normal(0, 20), "b"),
    prior(normal(0, 10), "sd"),
    prior(normal(0, 10), "Intercept", dpar = "sigma"),
    prior(normal(0, 10), dpar = "sigma")
)

# Running the model for Mean Theta
print("Running model for Mean Theta")
form <- brmsformula(Mean_Theta ~ 1 + Taxon + (1|Species) + (1|Population),
                    sigma ~ Taxon,
                    sparse = TRUE)
mod_mean_theta <- 
    brm_multiple(formula    = form,
                 data       = dfs_meta,
                 save_ranef = FALSE,
                 chains     = 2,
                 iter       = 3000,
                 warmup     = 2000,
                 thin       = 50, 
                 prior      = prior_n)
saveRDS(mod_mean_theta, file = "../mod_mean_theta.rds", version = 2)
rm(mod_mean_theta)
gc()

# Running the model for Mean Beta
print("Running model for Mean Beta")
form <- brmsformula(Mean_Beta ~ 1 + Taxon + (1|Species) + (1|Population),
                    sigma ~ Taxon,
                    sparse = TRUE)
mod_mean_beta <- 
    brm_multiple(formula    = form,
                 data       = dfs_meta,
                 save_ranef = FALSE,
                 chains     = 2,
                 iter       = 3000,
                 warmup     = 2000,
                 thin       = 50, 
                 prior      = prior_n)
saveRDS(mod_mean_beta, file = "../mod_mean_beta.rds", version = 2)
rm(mod_mean_beta)
gc()

# Running the model for SD Theta
print("Running model for SD Theta")
form <- brmsformula(SD_Theta  ~ 1 + Taxon + (1|Species) + (1|Population),
                    sigma ~ Taxon,
                    sparse = TRUE)
mod_sd_theta <- 
    brm_multiple(formula = form,
                 data   = dfs_meta,
                 save_ranef = FALSE,
                 chains = 2,
                 iter   = 3000,
                 warmup = 2000,
                 thin   = 50, 
                 prior  = prior_n)
saveRDS(mod_sd_theta, file = "../mod_sd_theta.rds", version = 2)
rm(mod_sd_theta)
gc()

# Running the model for SD Beta
print("Running model for SD Beta")
form <- brmsformula(SD_Beta ~ 1 + Taxon + (1|Species) + (1|Population),
                    sigma ~ Taxon,
                    sparse = TRUE)
mod_sd_beta <- 
    brm_multiple(formula = form,
                 data   = dfs_meta,
                 save_ranef = FALSE,
                 chains = 2,
                 iter   = 3000,
                 warmup = 2000,
                 thin   = 50, 
                 prior  = prior_n)
saveRDS(mod_sd_beta, file = "../mod_sd_beta.rds", version = 2)
rm(mod_sd_beta)
gc()

# Running the model for Omega
print("Running model for Omega")
form <- brmsformula(Omega ~ 1 + Taxon + (1|Species) + (1|Population),
                    sigma ~ Taxon,
                    sparse = TRUE)
mod_omega <- 
    brm_multiple(formula    = form,
                 data       = dfs_meta,
                 save_ranef = FALSE,
                 chains     = 2,
                 iter       = 3000,
                 warmup     = 2000,
                 thin       = 50, 
                 prior      = prior_n)
saveRDS(mod_omega, file = "../mod_omega.rds", version = 2)
rm(mod_omega)
gc()

## Formatting the five models and extracting outputs

est_meta <- 
    tibble(Param = c("Mean_Theta", "Mean_Beta", "SD_Theta", "SD_Beta", "Omega")) %>%
    mutate(Model = map(Param, ~ readRDS(str_glue("../mod_{str_to_lower(.)}.rds"))),
           Point_Est = map(Model, summarise_meta_brms)) %>%
    select(-Model) %>%
    unnest_legacy(Point_Est) %>%
    mutate(Taxon        = str_replace(Parameter, "mu_", ""),
           Parameter    = NULL) %>%
    select(Param, Taxon, everything())
saveRDS(est_meta, file = "../est_meta.rds", version = 2)

out_meta <-
    est_meta %>%
    mutate(Est = str_glue("{signif(Median, 3)} [{signif(Low, 2)}, {signif(Up, 2)}]")) %>%
    select(Param, Taxon, Est)
write_csv(out_meta, path = "../estimates_meta.csv")

## ---------------------------------- Plotting the graphs of the gradients

print("Generating the graphs")

# Formatting everything for use in ggplot2
graph_grad <-
    point_grad %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(Model = if_else(str_detect(Model_Year, "EXP"), "DIR.", "OPT.") %>%
                   as_factor()) %>%
    select(ID, Model, Mean_Beta, SD_Beta, SD_Beta_Nopl) %>%
    gather("Parameter", "Point_Estimates", -ID, -Model) %>% 
    mutate(Parameter = recode(Parameter,
                              Mean_Beta     = "E(beta)",
                              Var_Beta      = "V(beta)",
                              CV_Beta       = "CV(beta)",
                              SD_Beta       = "sigma(beta)",
                              SD_Beta_Nopl  = "sigma(beta['No pl.'])"),
           Point_Estimates = map(Point_Estimates, as_tibble),
           ID = as.factor(ID) %>% fct_rev()) %>% 
    unnest(Point_Estimates)

ylim <-
    graph_grad %>% 
    filter(str_detect(Parameter, "sigma")) %$%
    c(Low, Up) %>%
    range(na.rm = TRUE)

# Now creating the graph
p_mean <- 
    ggplot(graph_grad %>% filter(Parameter == "E(beta)")) +
    geom_linerange(aes(x = ID, ymin = Low, ymax = Up), size = 1) +
    geom_hline(yintercept = 0, colour = "grey", alpha = 0.5) +
    geom_point(aes(x = ID, y = Mode), colour = "#aa0000", size = 2) +
    geom_point(aes(x = ID, y = Median), colour = "#0055ff", size = 2) +
    ylab(bquote(E(beta))) + xlab("Dataset") +
    scale_colour_manual(values = c("#005500", "#55aa00")) +
    coord_flip() + 
    facet_wrap(~ Model, ncol = 1) +
    theme(text         = element_text(family = "Linux Biolinum O", size = 22),
          axis.text.y  = element_text(size = 8, family = "Linux Biolinum O"),
          axis.text.x  = element_text(size = 12, family = "Linux Biolinum O"),
          axis.title   = element_text(size = 22, family = "Linux Biolinum O"),
          legend.text  = element_text(size = 14),
          legend.title  = element_text(size = 16),
          strip.text   = element_text(family = "Linux Libertine O", size = 22),
          legend.position = "top")

p_sd <- 
    ggplot(graph_grad %>% filter(Parameter == "sigma(beta)")) +
    geom_linerange(aes(x = ID, ymin = Low, ymax = Up), size = 1) +
    geom_point(aes(x = ID, y = Mode), colour = "#aa0000", size = 2) +
    geom_point(aes(x = ID, y = Median), colour = "#0055ff", size = 2) +
    scale_y_log10(limits = ylim) +
    ylab(bquote(sigma(beta))) + xlab("Dataset") +
    scale_colour_manual(values = c("#005500", "#55aa00")) +
    coord_flip() + 
    facet_wrap(~ Model, ncol = 1) +
    theme(text         = element_text(family = "Linux Biolinum O", size = 22),
          axis.text.y  = element_text(size = 8, family = "Linux Biolinum O"),
          axis.text.x  = element_text(size = 12, family = "Linux Biolinum O"),
          axis.title   = element_text(size = 22, family = "Linux Biolinum O"),
          strip.text   = element_text(family = "Linux Libertine O", size = 22),
          legend.position = "none")

p_sd_nopl <- 
    ggplot(graph_grad %>% filter(Parameter == "sigma(beta['No pl.'])")) +
    geom_linerange(aes(x = ID, ymin = Low, ymax = Up), size = 1) +
    geom_point(aes(x = ID, y = Mode), colour = "#aa0000", size = 2) +
    geom_point(aes(x = ID, y = Median), colour = "#0055ff", size = 2) +
    scale_y_log10(limits = ylim) +
    ylab(bquote(sigma(beta['No pl.']))) + xlab("Dataset") +
    scale_colour_manual(values = c("#005500", "#55aa00")) +
    coord_flip() + 
    theme(text         = element_text(family = "Linux Biolinum O", size = 22),
          axis.text.y  = element_text(size = 8),
          axis.text.x  = element_text(size = 12),
          axis.title   = element_text(size = 22),
          strip.text   = element_text(family = "Linux Libertine O", size = 22),
          legend.position = "none")

# Now saving the graphs
cairo_pdf("../../Figures/Gradients.pdf", height = 10, width = 10)
p_mean + p_sd
dev.off()

# Now saving the graphs
cairo_pdf("../../Figures/Gradients_nopl.pdf", height = 15, width = 5)
p_sd / p_sd_nopl + plot_layout(heights = c(2,1))
dev.off()

# Changing the graphs for a talk
p_mean_talk <-
    p_mean +
    theme(text = element_text(family = "Noto Sans", size = 16),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
p_sd_talk <-
    p_sd +
    theme(text = element_text(family = "Noto Sans", size = 16),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
cairo_pdf("../../Figures/Gradients_preso.pdf", height = 10, width = 7)
p_mean_talk + p_sd_talk
dev.off()

## ------------------------------ Plotting the posterior distribution of the gradients

# Function to format betas
format_betas <-
    . %>%
    {suppressMessages(as_tibble(., .name_repair = "unique"))} %>%
    gather("Year", "Beta") %>%
    mutate(Year = str_remove(Year, "\\.\\.\\."))

# Formatting everything
graph_post <-
    all_grad %>%
    mutate(ID = str_c(DataID,
                      str_replace(Species, "^([A-Z])[a-z]+ ([a-z])[a-z ]+$", "\\1\\2"),
                      str_sub(Population, end = 3),
                      sep = "_"),
           Beta = future_map(Beta, format_betas, .progress = TRUE)) %>%
    select(ID, Beta) %>%
    unnest()

# Now creating the graph
p_post <-
    ggplot(graph_post) +
    geom_violin(aes(x = Year, y = Beta)) +
    facet_wrap(~ ID, scale = "free") +
    coord_flip() +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 10),
          axis.text.x  = element_text(size = 16),
          axis.title   = element_text(size = 18),
          strip.text   = element_text(family = "Linux Libertine O", size = 26))

# Now saving the graphs
cairo_pdf("../../Figures/Gradients_post.pdf", height = 25, width = 25)
plot(p_post)
dev.off()

## ------------------------------ Plotting prob. of sign variation

graph_sign <-
    point_grad %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(Model = if_else(str_detect(Model_Year, "EXP"), "DIR.", "OPT.") %>%
                   as_factor(),
           ID = fct_rev(ID)) %>%
    filter(Model == "OPT.") %>%
    select(ID, Var_Sign, Viable) %>% 
    unnest_wider(Var_Sign)

p_sign <- 
    ggplot(graph_sign) +
    geom_linerange(aes(x = ID, ymin = Low, ymax = Up), size = 1) +
    geom_hline(yintercept = 0, colour = "grey", alpha = 0.5) +
    geom_point(aes(x = ID, y = Median), colour = "#0055ff", size = 2) +
    ylab("Prop. of sign change") + xlab("Dataset") +
    scale_colour_manual(values = c("#005500", "#55aa00")) +
    coord_flip() + 
    theme(text         = element_text(family = "Linux Biolinum O", size = 22),
          axis.text.y  = element_text(size = 8, family = "Linux Biolinum O"),
          axis.text.x  = element_text(size = 12, family = "Linux Biolinum O"),
          axis.title   = element_text(size = 22, family = "Linux Biolinum O"),
          legend.text  = element_text(size = 14),
          legend.title  = element_text(size = 16),
          strip.text   = element_text(family = "Linux Libertine O", size = 22),
          legend.position = "top")

cairo_pdf("../../Figures/Gradients_VarSign.pdf", height = 6, width = 6)
plot(p_sign)
dev.off()

varsign <-
    point_grad %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(Model = if_else(str_detect(Model_Year, "EXP"), "DIR.", "OPT.") %>%
                   as_factor(),
           ID = fct_rev(ID)) %>%
    filter(Model == "OPT.") %>%
    transmute(ID = ID,
              Viable = Viable,
              Var_Sign  = map_dbl(Var_Sign, "Median"),
              Mean_Beta = map_dbl(Mean_Beta, "Median"),
              SD_Beta   = map_dbl(SD_Beta, "Median"),
              CV_Beta   = map_dbl(CV_Beta, "Median"))

varsign %$% cor(Var_Sign, Mean_Beta)
varsign %$% cor(Var_Sign, SD_Beta)
varsign %$% cor(Var_Sign, CV_Beta)

## ------------------------------ Plotting gradients and theta

others <- readRDS("../estimates.rds") %>%
          filter(Parameter %in% c("theta", "t_sigma", "omega")) %>%
          nest_legacy(Median, Mode, SE, Low, Up) %>%
          spread(Parameter, data) %>%
          select(everything(), Theta = theta, SD_Theta = t_sigma, Omega = omega)

all <- full_join(point_grad %>% select(-SD_Theta), others)

## Theta/Beta OPT
# Formatting "meta" estimates
formated_meta <-
    est_meta %>%
    filter(Param != "Omega") %>%
    mutate(Parameter = str_replace(Param, "Mean_Theta", "Theta"),
           Param     = NULL,
           Viable    = TRUE) %>%
    rename(ID = Taxon) %>%
    select(ID, Viable, Parameter, everything())

# Formatting Beta and Theta estimates, then adding "meta-analysis" estimates
all_graph_opt <-
    all %>%
    filter(Model_Year == "AR1") %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(ID = fct_rev(ID)) %>%
    select(ID, Theta, SD_Theta, Mean_Beta, SD_Beta, Viable) %>%
    mutate_at(vars(contains("Beta")), ~ map(., as_tibble)) %>%
    gather("Parameter", "Estimates", -ID, -Viable) %>%
    unnest(Estimates) %>%
    bind_rows(formated_meta) %>%
    mutate(ID = factor(ID, levels = c(levels(ids[["ID"]]), "Bird", "Mammal")) %>%
                fct_rev())

all_graph_opt <-
    all_graph_opt %>%
    mutate(
        Parameter_recode = Parameter %>%
                           recode_factor(
                               Theta = "E(θ)",
                               SD_Theta = "σ(θ)",
                               Mean_Beta = "E(β)",
                               SD_Beta = "σ(β)"
                           )
    )

p_opt <- 
    ggplot(all_graph_opt) +
    geom_linerange(aes(x = ID, ymin = Low, ymax = Up, colour = Viable), size = 1) +
    geom_point(aes(x = ID, y = Median, fill = Viable), shape = 21, size = 2) +
    coord_flip() +
    facet_wrap(~ Parameter_recode, scales = "free_x") +
    scale_colour_manual(values = c("grey70", "black"), guide = "none") +
    scale_fill_manual(values = c("#6d7494", "#0055ff"), guide = "none") +
    ylab("Estimates") + xlab("Data ID") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 8),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(family = "Linux Libertine O", size = 22))

cairo_pdf("../../Figures/Gradients_Theta_OPT.pdf", height = 10, width = 10)
plot(p_opt)
dev.off()

p_opt_talk <-
    ggplot(all_graph_opt %>% 
        filter(str_detect(Parameter, "Theta")) %>%
        filter(Viable)) +
    geom_linerange(aes(x = ID, ymin = Low, ymax = Up), size = 1) +
    geom_point(aes(x = ID, y = Median), size = 2) +
    coord_flip() +
    facet_wrap(~ Parameter_recode, scales = "free_x", ncol = 1) +
    ylab("Estimates") + xlab("Data ID") +
    theme(text = element_text(family = "Noto Sans", size = 16),
          strip.text   = element_text(family = "Noto Sans", size = 16),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())

cairo_pdf("../../Figures/Theta_preso.pdf", height = 10, width = 5)
plot(p_opt_talk)
dev.off()

## Theta/Beta OPT with Weigths
weights <-
    loo %>%
    left_join(ids) %>%
    select(ID, Model_Year, LOO_Weight) %>%
    group_by(ID) %>%
    summarise(Weight_Fluct = sum(LOO_Weight[str_detect(Model_Year, "(AR1)|(IID)")]),
              Weight_Opt = sum(LOO_Weight[Model_Year %in% c("AR1", "IID", "NULL")]))


all_graph_weights <-
    all %>%
    filter(Model_Year == "AR1") %>%
    select(-ID) %>%
    left_join(ids) %>%
    left_join(weights) %>%
    mutate(ID = fct_rev(ID)) %>%
    select(ID, Theta, SD_Theta, Mean_Beta, SD_Beta, Weight_Fluct, Weight_Opt) %>%
    mutate_at(vars(contains("Beta")), ~ map(., as_tibble)) %>%
    gather("Parameter", "Estimates", -ID, -Weight_Fluct, -Weight_Opt) %>%
    unnest(Estimates) %>%
    bind_rows(formated_meta) %>%
    mutate(ID = factor(ID, levels = c(levels(ids[["ID"]]), "Bird", "Mammal")) %>%
                fct_rev())

all_graph_weights <-
    all_graph_weights %>%
    mutate(
        Parameter_recode = Parameter %>%
                           recode_factor(
                               Theta = "E(θ)",
                               SD_Theta = "σ[θ]",
                               Mean_Beta = "E(β)",
                               SD_Beta = "σ[β]"
                           )
    )

p_weights <- 
    ggplot(all_graph_weights) +
    geom_hline(yintercept = 0, colour = "grey50") + 
    geom_linerange(data = all_graph_weights %>% filter(Parameter %in% c("SD_Theta", "SD_Beta")),
                   aes(x        = ID,
                       ymin     = Low,
                       ymax     = Up,
                       colour   = Weight_Fluct,
                       size     = is.na(Weight_Opt))) +
    geom_point(data = all_graph_weights %>% filter(Parameter %in% c("SD_Theta", "SD_Beta")),
               aes(x = ID, y = Median, colour = Weight_Fluct), size = 2) +
    scale_colour_gradient(low = "#000000",
                          high = "#ff8000",
                          name = "Fluctuation support",
                          na.value = "#00aa00",
                          guide = guide_colourbar(order = 2),
                          limits = c(0,1)) +
    new_scale_colour() +
    geom_linerange(data = all_graph_weights %>% filter(Parameter %in% c("Theta")),
                   aes(x        = ID,
                       ymin     = Low,
                       ymax     = Up,
                       colour   = Weight_Opt,
                       size     = is.na(Weight_Opt))) +
    geom_point(data = all_graph_weights %>% filter(Parameter %in% c("Theta")),
               aes(x = ID, y = Median, colour = Weight_Opt), size = 2) +
    geom_linerange(data = all_graph_weights %>% filter(Parameter %in% c("Mean_Beta")),
                   aes(x        = ID,
                       ymin     = Low,
                       ymax     = Up,
                       colour   = Weight_Opt,
                       size     = is.na(Weight_Opt))) +
    geom_point(data = all_graph_weights %>% filter(Parameter %in% c("Mean_Beta")),
               aes(x = ID, y = Median, colour = Weight_Opt), size = 2) +
    scale_colour_gradient(low = "#000000",
                          high = "#0080ff",
                          name = "Optimum support",
                          guide = guide_colourbar(order = 1),
                          na.value = "#00aa00",
                          limits = c(0,1)) +
    geom_point(data = all_graph_weights %>% filter(ID %in% c("Bird", "Mammal")),
               aes(x = ID, y = Median),
               colour = "#00aa00",
               size = 3,
               shape = 22,
               fill = "#00aa00") +
    scale_size_discrete(guide = "none", range = c(0.8, 1.8)) +
    scale_x_discrete(expand = expansion(add = 1)) +
    coord_flip() +
    facet_wrap(~ Parameter_recode, scales = "free_x", labeller = label_parsed) +
    ylab("Estimates") + xlab("Data ID") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 10),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 20),
          legend.title  = element_text(size = 20, vjust = 1),
          legend.text   = element_text(size = 14, angle = 45, vjust = 1.5, hjust = 1),
          legend.box.spacing = element_blank(),
          strip.text   = element_text(family = "Linux Libertine O", size = 22),
          legend.position = "top")

cairo_pdf("../../Figures/Gradients_Theta_OPT_weights.pdf", height = 11, width = 10)
plot(p_weights)
dev.off()

## Theta/Beta on the same graph
all_graph_same <-
    all %>%
    filter(Model_Year == "AR1", Viable) %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(ID = fct_rev(ID)) %>%
    select(ID, Theta, SD_Theta, Omega, Mean_Beta, SD_Beta, Viable) %>%
    mutate_at(vars(contains("Beta")), ~ map(., as_tibble)) %>%
    gather("Parameter", "Estimates", -ID, -Viable) %>%
    unnest(Estimates) %>%
    filter(Parameter %in% c("Theta", "Mean_Beta", "Omega")) %>%
    select(-Mode, -SE) %>%
    pivot_wider(names_from = Parameter, values_from = c(Median, Low, Up))

p_same <-
    ggplot(all_graph_same,
           aes(x = Median_Theta, y = Median_Mean_Beta, label = ID, colour = Median_Omega)) +
    geom_point() +
    geom_text_repel(family = "Linux Biolinum O", box.padding = 0.3) +
    xlab("E(θ)") + ylab("E(β)") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(family = "Linux Libertine O", size = 8),
          axis.text.x  = element_text(family = "Linux Libertine O", size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(family = "Linux Biolinum O", size = 22)) +
    scale_colour_gradient(low = "#aa0000", high = "#ffaa00", name = "ω")

cairo_pdf("../../Figures/Gradients_Theta_Sameplot.pdf", height = 6, width = 7)
plot(p_same)
dev.off()

## Theta/Beta DIR
all_graph_dir <-
    all %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(ID = fct_rev(ID)) %>%
    select(ID, Model_Year, Theta, SD_Theta, Mean_Beta, SD_Beta, Viable) %>%
    mutate_at(vars(contains("Beta")), ~ map(., as_tibble)) %>%
    gather("Parameter", "Estimates", -ID, -Viable, -Model_Year) %>%
    unnest(Estimates) %>%
    filter(
        (str_detect(Parameter, "Theta") & Model_Year == "AR1") |
            (str_detect(Parameter, "Beta") & Model_Year == "EXPAR1")
    ) %>%
    mutate(Viable = if_else(Model_Year == "EXPAR1", TRUE, Viable))

all_graph_dir <-
    all_graph_dir %>%
    mutate(
        Parameter = Parameter %>%
                    recode_factor(
                        Theta = "E(θ)",
                        SD_Theta = "σ(θ)",
                        Mean_Beta = "E(β)",
                        SD_Beta = "σ(β)"
                    )
    )

p_dir <- 
    ggplot(all_graph_dir) +
    geom_linerange(aes(x = ID, ymin = Low, ymax = Up, colour = Viable), size = 1) +
    geom_point(aes(x = ID, y = Median, fill = Viable), shape = 21, size = 2) +
    coord_flip() +
    facet_wrap(~ Parameter, scales = "free_x") +
    scale_colour_manual(values = c("grey70", "black"), guide = "none") +
    scale_fill_manual(values = c("#6d7494", "#0055ff"), guide = "none") +
    ylab("Estimates") + xlab("Data ID") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 8),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(family = "Linux Libertine O", size = 22))

cairo_pdf("../../Figures/Gradients_Theta_DIR.pdf", height = 10, width = 10)
plot(p_dir)
dev.off()

## Compare DIR/OPT

graph_comp <-
    point_grad %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(ID = fct_rev(ID)) %>%
    select(ID, Model_Year, Mean_Beta, SD_Beta) %>%
    mutate_at(vars(contains("Beta")), ~ map_dbl(., "Median")) %>%
    gather("Parameter", "Value", -ID, -Model_Year) %>%
    mutate(Parameter = str_remove(Parameter, "_Beta")) %>%
    spread(Model_Year, Value)

graph_comp %>% filter(Parameter == "Mean") %$% cor(AR1, EXPAR1) # 0.98
lm_mean <-
    graph_comp %>%
    filter(Parameter == "Mean") %$%
    lm(AR1 ~ EXPAR1) %>%
    summary() %>%
    {c(.[["coefficients"]]["EXPAR1", "Estimate"], .[["r.squared"]])} %>%
    round(digits = 2)
# slope = 0.92, R² = 0.96

graph_comp %>% filter(Parameter == "SD") %$% cor(AR1, EXPAR1)   # 0.89
lm_sd <-
    graph_comp %>%
    filter(Parameter == "SD") %$%
    lm(AR1 ~ EXPAR1) %>%
    summary() %>%
    {c(.[["coefficients"]]["EXPAR1", "Estimate"], .[["r.squared"]])} %>%
    round(digits = 2) 
# slope = 0.98, R² = 0.79

df_annot <-
    tibble(Parameter = c("Mean", "SD"),
           Slope     = c(lm_mean[1], lm_sd[1]),
           R2        = c(lm_mean[2], lm_sd[2])) %>%
    mutate(Text = str_glue("b = {Slope},\n R² = {R2}"),
           X    = c(-0.52, 0.18),
           Y    = c(-0.37, 0.15))

p_comp <-
    ggplot(graph_comp) +
    geom_abline(intercept = 0, slope = 1, colour = "grey50") +
    geom_smooth(aes(x = EXPAR1, y = AR1),
                method  = "lm",
                colour  = "#0055ff",
                se      = FALSE) +
    geom_point(aes(x = EXPAR1, y = AR1)) +
    geom_text_repel(aes(x = EXPAR1, y = AR1, label = ID),
                    segment.alpha = 0.5) +
    geom_text(data  = as.data.frame(df_annot),
              aes(label = Text, x = X, y = Y),
              hjust     = 0,
              colour    = "#0055ff",
              family    = "Linux Libertine O",
              lineheight= 0.8,
              size      = 5) +
    xlab("β from FluctCorrDir") + ylab("β from FluctCorrOpt") +
    facet_wrap(~ Parameter, scale = "free") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(family = "Linux Biolinum O", size = 22))

cairo_pdf("../../Figures/Gradients_DIR_VS_OPT.pdf", height = 6, width = 12)
plot(p_comp)
dev.off()

## SD Theta - Zbar

graph_sdtz <-
    opt %>%
    mutate(Zbar     = map2(Zbar, Theta, ~ {matrix(.x, ncol = length(.x), nrow = nrow(.y), byrow = TRUE)}),
           SD_Zbar  = map_dbl(Zbar, ~ sd(., na.rm = TRUE)),
           CorrTZ   = future_map2(Theta, Zbar,
                                  ~ apply(array(c(.x, .y), dim = c(dim(.x), 2)),
                                          1,
                                          function(mat) {cor(mat[ , 1], mat[ , 2], use = "complete.obs")}))) %>%
    select(-SD_Theta, -Omega) %>%
    left_join(others %>% select(DataID, Species, Population, SD_Theta, Omega)) %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(ID = fct_rev(ID)) %>%
    select(ID, SD_Theta, SD_Zbar, CorrTZ) %>%
    mutate(CorrTZ   = map_dbl(CorrTZ, ~ median(., na.rm = TRUE)),
           SD_Theta = map_dbl(SD_Theta, "Median"))

p_sdtz <-
    ggplot(graph_sdtz, aes(x = SD_Zbar, y = SD_Theta, colour = CorrTZ, label = ID)) +
    geom_point() +
    geom_text_repel(family = "Linux Biolinum O",
                    box.padding = 0.3,
                    min.segment.length = 0.35) +
    scale_colour_gradient2(low = "#0000ff", high = "#ff0000", mid = "black", midpoint = 0,
                           name = bquote("ρ[θ,z]")) +
    xlab(bquote(paste("σ(", bar(z) ,")"))) + ylab("σ(θ)") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(family = "Linux Libertine O", size = 22))

cairo_pdf("../../Figures/Gradients_SD_Theta_Z.pdf", height = 6, width = 6)
plot(p_sdtz)
dev.off()

## SD with or without plasticity

graph_plast <-
    point_grad %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(ID = fct_rev(ID)) %>%
    filter(Model_Year == "AR1", Viable) %>%
    select(ID, SD_Beta, SD_Beta_Nopl, SD_Zbar, SD_Theta) %>%
    mutate(SD_Beta = map_dbl(SD_Beta, "Median"),
           SD_Beta_Nopl = map_dbl(SD_Beta_Nopl, "Median"),
           SD_Theta     = map_dbl(SD_Theta, median, na.rm = TRUE),
           Ratio        = SD_Zbar / SD_Theta)

p_plast <-
    ggplot(graph_plast, aes(x = SD_Beta_Nopl, y = SD_Beta, colour = Ratio, label = ID)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey60", size = 1) +
    geom_point() +
    geom_text_repel(family = "Linux Biolinum O",
                    box.padding = 0.3,
                    min.segment.length = 0.35) +
    xlab(bquote(paste(frac(σ[θ], ω^2 + 1)))) + ylab(bquote(paste("σ(",β[t],")"))) + 
    scale_x_log10() + scale_y_log10() +
    scale_colour_gradient(low = "#000000", high = "#ff0000",
                          name = bquote(frac(paste("σ(", bar(z) ,")"), σ[θ]))) +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(family = "Linux Libertine O", size = 22))

cairo_pdf("../../Figures/Gradients_Plasticity.pdf", height = 6, width = 7)
plot(p_plast)
dev.off()

## SD with or without tracking

graph_track <-
    point_grad %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(ID = fct_rev(ID)) %>%
    filter(Model_Year == "AR1", Viable) %>%
    select(ID, SD_Beta, SD_Beta_Notr, SD_Theta, SD_Zbar, CorrTZ) %>%
    mutate(SD_Beta = map_dbl(SD_Beta, "Median"),
           SD_Beta_Notr = map_dbl(SD_Beta_Notr, "Median"),
           CorrTZ       = map_dbl(CorrTZ, median, na.rm = TRUE),
           SD_Theta     = map_dbl(SD_Theta, median, na.rm = TRUE),
           Ratio        = SD_Zbar / SD_Theta)

p_track <-
    ggplot(graph_track, aes(x = SD_Beta_Notr, y = SD_Beta, colour = CorrTZ, label = ID)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey60", size = 1) +
    geom_point() +
    geom_text_repel(family = "Linux Biolinum O",
                    box.padding = 0.3,
                    min.segment.length = 0.35) +
    xlab(bquote(paste(frac(sqrt(σ[θ]^2 + paste("σ(", bar(z) ,")")^2), ω^2 + 1)))) + ylab(bquote(paste("σ(",β[t],")"))) + 
    scale_x_log10() + scale_y_log10() +
    scale_colour_gradient2(low = "#0000ff", high = "#ff0000", mid = "black", midpoint = 0,
                           name = bquote("ρ[θ,z]")) +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(family = "Linux Libertine O", size = 22))

cairo_pdf("../../Figures/Gradients_Tracking.pdf", height = 6, width = 7)
plot(p_track)
dev.off()

## SD complete and final graph

graph_all <-
    point_grad %>%
    select(-ID) %>%
    left_join(ids) %>%
    mutate(ID = fct_rev(ID)) %>%
    filter(Model_Year == "AR1", Viable) %>%
    select(ID, SD_Beta, SD_Beta_Nofl, SD_Beta_Nopl, SD_Beta_Notr, CorrTZ, SD_Zbar, SD_Theta) %>%
    mutate(SD_Beta      = map_dbl(SD_Beta, "Median"),
           SD_Beta_Nofl = map_dbl(SD_Beta_Nofl, "Median"),
           SD_Beta_Nopl = map_dbl(SD_Beta_Nopl, "Median"),
           SD_Beta_Notr = map_dbl(SD_Beta_Notr, "Median"),
           CorrTZ       = map_dbl(CorrTZ, median, na.rm = TRUE),
           SD_Theta     = map_dbl(SD_Theta, median, na.rm = TRUE),
           Ratio        = SD_Zbar / SD_Theta,
           Nudge        = if_else(CorrTZ > 0.25, -0.25, -0.05),
           Nudge        = if_else(ID %in% c("Cca1", "Pdo", "Cca6", "Cci1", "Pma1", "Cel", "Rta", "Oca", "Cca8"), 0.2, Nudge),
           Nudge        = if_else(ID %in% c("Pma1", "Cci2"), -0.15, Nudge))

p_all <-
    ggplot(graph_all) +
    geom_abline(intercept = 0, slope = 1, colour = "grey60", size = 1) +
    geom_text_repel(aes(x = SD_Beta_Nopl, y = SD_Beta, label = ID),
                    colour = "grey50",
                    segment.colour = "grey50",
                    segment.alpha  = 0.35,
                    family = "Linux Biolinum O",
                    box.padding = 0.3,
                    min.segment.length = 0.35,
                    nudge_y = graph_all[["Nudge"]]) +
    geom_point(aes(x = SD_Beta_Nopl, y = SD_Beta_Notr), shape = 4) +
    geom_point(aes(x = SD_Beta_Nopl, y = SD_Beta)) +
    geom_segment(aes(x      = SD_Beta_Nopl,
                     y      = SD_Beta_Notr,
                     xend   = SD_Beta_Nopl,
                     yend   = 1.03 * SD_Beta,
                     colour = CorrTZ^2),
                 arrow      = arrow(angle = 25, type = "closed", length = unit(0.011, "npc")),
                 linejoin   = 'round') +
    xlab(bquote(σ[paste("β (only optimum fluctuates)")])) +
    ylab(bquote(σ[paste("β (optimum and mean trait fluctuate)")])) + 
    annotate(geom = "rect", xmin = 0.003, xmax = 0.012, ymin = 0.075, ymax = 0.48,
             fill = "white", colour = "black", linejoin = "round") +
    annotate(geom = "point", x = 0.006, y = 0.25, shape = 4) +
    annotate(geom = "point", x = 0.006, y = 0.15) +
    annotate(geom = "segment", x = 0.006, xend = 0.006, y = 0.25, yend = 1.03 * 0.15,
             arrow      = arrow(angle = 25, type = "closed", length = unit(0.011, "npc")),
             linejoin   = 'round', colour = "#cf5c42") +
    annotate(geom = "text", x = 0.006, y = 0.29,
             label = "Estimate\nwithout\ntracking", vjust = 0, hjust = 0.5,
             lineheight = 0.7, family = "Linux Biolinum O", size = 5) +
    annotate(geom = "text", x = 0.006, y = 0.135,
             label = "Actual\nestimate\n(incl. tracking)", vjust = 1, hjust = 0.5,
             lineheight = 0.7, family = "Linux Biolinum O", size = 5) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_gradient(low = "grey75", high = "#aa0000",
                           name = bquote(ρ[paste(bar(z),",θ")]^2)) +
    theme(text          = element_text(family = "Linux Biolinum O"),
          axis.text.y   = element_text(size = 18),
          axis.text.x   = element_text(size = 18),
          axis.title    = element_text(size = 26),
          legend.position = "top",
          legend.title  = element_text(size = 20),
          legend.text   = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
          strip.text    = element_text(family = "Linux Libertine O", size = 22))

cairo_pdf("../../Figures/Gradients_Final.pdf", height = 7, width = 6)
plot(p_all)
dev.off()

## Some R²
lm(SD_Beta ~ SD_Beta_Nopl, data = graph_all) %>% summary() # R² = 0.98
lm(SD_Beta ~ SD_Beta_Nofl, data = graph_all) %>% summary() # R² = 0.11
