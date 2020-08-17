library(tidyverse)

filter  <- dplyr::filter
lag     <- dplyr::lag
extract <- rstan::extract

## ---------------- Loading the data

# List of IDs
ids  <- readRDS("../FinalIDs.rds")

# Estimates
est <- 
    readRDS("../estimates.rds") %>%
    select(-Viable)

# Gradients
grad <- 
    readRDS("../gradients.rds") %>%
    select(DataID, Species, Population, Dist, Model_Year, Mean_Beta, SD_Beta, Keep_Theta) %>%
    pivot_longer(contains("Beta"), names_to = "Parameter", values_to = "Est") %>%
    unnest_wider(Est) %>%
    mutate(Parameter = recode(Parameter, Mean_Beta = "beta", SD_Beta = "sigma_beta"))

## ---------------- Saving a table with all the estimates

csv_tbl <-  bind_rows(est, grad)

# Formatting the data to save the published CSV
csv_tbl <-
    csv_tbl %>%
    left_join(ids) %>%
    select(ID, Dist, Model_Year, Parameter:Up, Keep_Theta) %>%
    mutate(Model        = recode(Model_Year,
                                 AR1 = "FluctCorrOpt",
                                 EXPAR1 = "FluctCorrDir"),
           Model_Year   = NULL,
           Parameter = factor(Parameter,
                              levels = c("theta", "theta_unstd", "t_sigma", "t_sigma_unstd",
                                         "log_Wmax", "logit_Wmax",
                                         "log_w_sigma", "logit_w_sigma",
                                         "slope", "sl_sigma", "intercept", "int_sigma",
                                         "omega", "phi", "p_zi", "ind_sigma",
                                         "beta", "sigma_beta")),
           Parameter    = recode(Parameter,
                                 t_sigma       = "sigma_theta",
                                 t_sigma_unstd = "sigma_theta_unstd",
                                 log_w_sigma   = "sigma_log_Wmax",
                                 intercept     = "a",
                                 int_sigma     = "sigma_a",
                                 slope         = "b",
                                 sl_sigma      = "sigma_b",
                                 ind_sigma     = "sigma_ind")) %>%
    rename(Post_median  = Median, 
           Post_mode    = Mode,
           Post_SE      = SE,
           CI95_Low     = Low,
           CI95_Up      = Up,
           Majority_Opt = Keep_Theta) %>%
    select(ID, Dist, Model, Parameter, everything()) %>%
    arrange(ID, Model, Parameter) %>%
    mutate_if(is.numeric, ~ signif(., digits = 3))

write_csv(csv_tbl, path = "../estimates.csv")
