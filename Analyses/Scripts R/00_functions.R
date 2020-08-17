library(tidyverse)
library(magrittr)

## --------------------------------- Solving possible conflicts
filter <- dplyr::filter
lag    <- dplyr::lag

## --------------------------------- Small functions

## Standardising function that is simpler than "scale"
simple_scale <-
    . %>%
    {(. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)}

## Minmax standardising function
scale_minmax <-
    . %>%
    {(. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))}

## Is a vector composed of a unique value
is_unique <-
    . %>%
    {length(unique(.)) == 1}

## Small functions to restart the progress bar
restart_pb <- function(total, force = TRUE, format = "[:bar] :percent ETA: :eta") {
    pb <<- progress::progress_bar$new(total   = total,
                                      force   = force,
                                      format  = format)
}

## Compute variations in sign
compute_var_sign <- function(vec) {
    sum(sign(vec) != sign(median(vec)), na.rm = TRUE) / length(which(!is.na(vec)))
}

## --------------------------------- Formatting and data wrangling functions

## Function to format data to analyse models
# Args: - tbl : tbl containing all the data (e.g. "alldata")
# Value: a nested tbl for each "dataset"
format_to_mods <- function(tbl) {
    
    # Getting the current year for transforming Pheno to a numerical
    current_year <- lubridate::year(lubridate::today())
    
    # Starting the formatting
    tbl <-
        tbl %>%    
        # Grouping by dataset
        group_by(DataID, Species, Population) %>%
        # Need rownames to create dummy IDs
        rownames_to_column() %>%
        # Creating a few variables (remember tbl is grouped by dataset)
        mutate(
            # Creating dummy IDs for missing female IDs
            ID   = if_else(is.na(ID),    
                           stringr::str_glue("Dummy_{rowname}") %>% as.character(),
                           ID),
            # Phenology to a numerical value
            Pheno_num   = {Pheno - lubridate::ymd(stringr::str_glue("{current_year}-01-01"))} %>%
                           as.double(),
            SD_Within   = lm(Pheno_num ~ as.factor(Year)) %>% 
                          sigma(),
            # Scaling phenology
            Pheno_scale = (Pheno_num - mean(Pheno_num, na.rm = TRUE)) / SD_Within,
            # Pheno squared
            Pheno_sq = Pheno_scale^2
        ) %>%
        # Remove rownames
        select(-rowname)
        
    # Now nesting and additional formatting
    tbl <-
        tbl %>%
        # Nesting the data within each "dataset"
        nest_legacy() %>%
        # Renaming the "data" column to avoid issue with a "data" variable
        rename(Data = data) %>%
        # Getting information for when to use binomial distribution
        mutate(Dist = Data %>% 
                      map_lgl(~ {max(.[["Fitness"]]) == 1}) %>%
                      if_else("binomial", "poisson")) %>%
        # Is ID missing
        mutate(Missing_ID = map_lgl(Data,
                                ~ all(str_detect(.[["ID"]], "Dummy")))) %>%
        # Creating a factor for years for each dataset (throw warnings in mutate above)
        mutate(Data = Data %>%
                      map(~ mutate(., Year_fac = factor(Year, levels = min(Year):max(Year))))) %>%
        # Organizing the tbl
        select(DataID, Species, Population, Dist, Data)
    
    return(tbl)
}

## Function to create formulas depending on some paramters
# Args: - quad: should the model include random effects (Yes/No)
#       - random: type of random slope (None/IID/AR1)
# Value: a formula to be passed on to INLA
generate_formulas_inla <- function(quad, random) {
    string <- 
        "Fitness ~ " %>%
        stringr::str_c(if_else(random == "None",
                      "",
                      " Year_fac +")) %>%
        stringr::str_c(if_else(quad == "Yes",
              " Pheno_scale + Pheno_sq",
              " Pheno_scale")) %>%
        stringr::str_c(case_when(
            random == "None"    ~ '',
            random == "IID"     ~ '+ f(Year_fac, Pheno_scale, model = "iid")',
            random == "AR1"     ~ '+ f(Year_fac, Pheno_scale, model = "ar1")'
        ))
    as.formula(string)
}

## Function to transform INLA parameters into STAN parameters
# Args: - model: An INLA model
# Value: a list of initial values to be passed to STAN
inla_to_init_stan <- function(model) {
        
    # Getting the estimates
    if (length(model[["summary.random"]]) == 0) {
        est <- tibble(Params = c(rownames(model[["summary.fixed"]])),
                      Values = c(model[["summary.fixed"]][ , "mean"]))
    } else {
        est <- tibble(Params = c(rownames(model[["summary.fixed"]]),
                                 model[["summary.random"]][[1]][ , "ID"]),
                      Values = c(model[["summary.fixed"]][ , "mean"],
                                 model[["summary.random"]][[1]][ , "mean"]))
    }
    
    # Formatting the whole mess
    est <-
        est %>%
        mutate(Parameter = Params %>%
                           # Starting with Year_fac are the intercept factors (hence b0)
                           stringr::str_replace("^Year_fac[0-9]+$", "B0") %>%
                           # Just years are the random slope estimates (hence b1)
                           stringr::str_replace("^[0-9]+$", "B1") %>%
                           # Renaming intercept (b0), 1st order slope (b1) and 2nd order slope (b2)
                           recode(`(Intercept)` = "B0",
                                  `Pheno_scale` = "B1",
                                  `Pheno_sq` = "B2"),
               # Getting the years when available
               Year      = Params %>%
                           stringr::str_extract("[0-9]+$") %>%
                           {if_else(is.na(.), "Overall",.)},
               Params    = NULL) %>%
        # "Overall" values should be common to everyone for it to be tidy
        mutate(Overall_B0 = Values[Parameter == "B0" & Year == "Overall"],
               Overall_B1 = Values[Parameter == "B1" & Year == "Overall"],
               B2         = Values[Parameter == "B2" & Year == "Overall"])
    
    # Need to have different strategies when there are no year-to-year variation
    if (est %$% all(Year == "Overall")) {
        # Simply use the available values and return them
        est <-
            est %>%
            select(Parameter, Values) %>%
            spread(Parameter, Values)
    } else {
        # In that case, we need to compute everything year by year
        est <-
            est %>%
            # Now we can remove the "Overall" values
            dplyr::filter(Year != "Overall") %>%
            # Let's take care of per-year b0 and b1 now
            spread(Parameter, Values) %>% 
            mutate(B0 = if_else(is.na(B0), 0, B0),  # First value was blended into the intercept by R
                   B0 = B0 + Overall_B0,    # Contrasts to the overall effect
                   B1 = B1 + Overall_B1,    # Contrasts to the overall effect
                   Overall_B0 = NULL,
                   Overall_B1 = NULL,
                   Year = as.integer(Year)) %>%
            # Reordering nicely
            select(Year, B0, B1, B2)
    }
    
    # Need to account for when B2 is positive
    est <-
        est %>%
        mutate(B2 = if_else(B2 > 0,
                            -0.001,
                            B2))
    
    # Formatting to explicit model
    est <- 
        est %>%
        mutate(Explicit = pmap(list(B0, B1, B2), ~ glmm_to_explicit(..1, ..2, ..3))) %>% 
        unnest_legacy() %>%
        select(-B0, -B1, -B2)
    
    # Obtaining p_zi if applicable
    detected <- stringr::str_detect(rownames(model[["summary.hyperpar"]]), "zero-probability")
    if (any(detected)) {
        ZI <- model[["summary.hyperpar"]][detected, "mean"]
    } else {
        ZI <- NULL
    }
    
    # Obtaining the auto-correlation when applicable
    detected <- stringr::str_detect(rownames(model[["summary.hyperpar"]]), "Rho")
    if (any(detected)) {
        Phi <- model[["summary.hyperpar"]][detected, "mean"]
    } else {
        Phi <- NULL
    }
    
    # Outputting the results as a list
    init <- list(omega  = mean(est[["Omega"]]),
                 p_zi   = ZI,
                 phi    = Phi,
                 log_Wmax   = log(mean(est[["Wmax"]])),
                 theta  = mean(est[["Theta"]]),
                 t_intercept    = mean(est[["Theta"]]),
                 t_year         = est[["Theta"]],
                 t_sigma        = sd(est[["Theta"]]),
                 log_w_intercept    = log(mean(est[["Wmax"]])),
                 log_w_year         = log(est[["Wmax"]]),
                 log_w_sigma        = log(sd(est[["Wmax"]])))
    
    if (is.na(init[["t_sigma"]])) { init[["t_sigma"]] <- NULL }
    if (is.na(init[["w_sigma"]])) { init[["w_sigma"]] <- NULL }
    
    return(init)
}

## Function to set up a vector of parameters to keep in some STAN output
# Args: - mod_year: the model used for between-year variations (NULL/IID/AR1)
#       - dist: distribution used for the data (e.g. is it ZIPoisson?)
# Value: a vector of parameters to keep on some STAN outputs
get_params_to_keep <- function(mod_year, dist) {   
    # Parameters depending on mod_year
    if (mod_year == "NULL") {
        out <- c("omega", "theta", "log_Wmax", "ind_sigma")
    } else if (mod_year == "VARINT") {
        out <- c("omega", "theta", "log_w_intercept", "log_w_sigma", "ind_sigma")
    } else if (mod_year == "EXPNULL") {
        out <- c("intercept", "slope", "ind_sigma")
    } else if (mod_year == "EXPVARINT") {
        out <- c("intercept", "int_sigma", "slope", "ind_sigma")
    } else if (mod_year == "EXPIID") {
        out <- c("intercept", "int_sigma", "slope", "sl_sigma", "ind_sigma")
    } else if (mod_year == "EXPAR1") {
        out <-c("intercept", "int_sigma", "slope", "sl_sigma", "phi", "ind_sigma")
    } else if (mod_year == "IID") {
        out <- c("omega", "t_intercept", "t_sigma", "log_w_intercept", "log_w_sigma", "ind_sigma")
    } else if (mod_year == "AR1") {
        out <- c("omega", "t_intercept", "t_sigma", "log_w_intercept", "log_w_sigma", "phi", "ind_sigma")
    } else if (mod_year == "AR1_withtrend") {
        out <- c("omega", "t_intercept", "t_sigma", "log_w_intercept", "log_w_sigma", "trend", "phi", "ind_sigma")
    } else if (mod_year == "FLAT") {
        out <- c("log_w_intercept", "log_w_sigma", "ind_sigma")
    }
    
    # Finally, adding the zero-inflation for ZIPoisson
    if (dist == "ZIPoisson") {
        out <- c(out, "p_zi")
    }
    
    return(out)
}

## Function to add values of Wmax on the log scale to a mcmc
# Args: - mcmc: a mcmc object
# Value: a mcmc object with log_w_intercept and log_w_sigma added to it
add_logWmax_to_mcmc <- function(mcmc) {
    if ("Wmax" %in% coda::varnames(mcmc)) {
        out <- cbind(mcmc,
                     log_Wmax = log(mcmc[ , "Wmax"]))
        out <- out[ , -grep("^Wmax$", colnames(out))]
        out <- coda::as.mcmc(out)
        return(out)
    } else if (coda::varnames(mcmc) %>% str_detect("^w_") %>% any()) {
        tbl <- transform_Wlog(mcmc[ , "w_intercept"], mcmc[ , "w_sigma"])
        out <- cbind(mcmc,
                     log_w_intercept = tbl[["log_w_intercept"]],
                     log_w_sigma     = tbl[["log_w_sigma"]])
        out <- out[ , -grep("^w_", colnames(out))]
        out <- coda::as.mcmc(out)
        return(out)
    } else {
        return(mcmc)
    }
}

## Function to obtain average lay date each year from a dataset
# Args: - data: a data tbl containing the data used to fit the model (a tbl)
#       - column: a value for the column to get the data from (a string)
# Value: a vector of zbar values for each year (including missing ones)
compute_zbar <- function(data, column) {
    
    # Getting the mean phenotype for each year and complete missing years with NA
    zbar <-
        data %>% 
        group_by(Year) %>% 
        summarise(Zbar = mean(.data[[column]])) %>%
        complete(Year = full_seq(Year, 1)) %$%  # NA are attributed to missing years
        {Zbar}
    
    return(zbar)
}

## Function to obtain differences in lay date for each individual between years from a dataset
# Args: - tbl: a data tbl containing the data used to fit the model (a tbl)
#       - column: a value for the column to get the data from (a string)
# Value: a vector of zbar values for each year (including missing ones)
compute_delta_zbar <- function(tbl, column) {
    
    # Getting the mean phenotype for each year and complete missing years with NA
    deltas <-
        tbl %>% 
        nest_legacy(-Year) %>%
        mutate(ID = map(data, ~ .[["ID"]]),
               # Saving the total number of individuals
               N_ID = map_int(ID, length),
               # Getting individuals in common between consecutive years
               CommonID = map2(ID, lead(ID), ~ .x[.x %in% .y]),
               # Saving the number of such individuals
               N_CommonID = map_int(CommonID, length),
               # Computing the Zbar for each on all individuals (for possible checks)
               Zbar = map_dbl(data,
                              ~ mean(.[[column]])),
               # Computing the Zbar on the "in common" individuals (for possible checks)
               CommonZbar = map2_dbl(data, CommonID,
                                     ~ {filter(.x, ID %in% .y) %$% mean(Pheno_scale)}),
               # Computing the differences in Zbar for all individuals (for possible checks)
               DeltaZbar = map2_dbl(data, lead(data),
                                    ~ {mean(pluck(.y, column, .default = NA)) - 
                                       mean(.x[[column]])}),
               # Computing the differences in Zbar for all individuals (focal result)
               DeltaCommonZbar = pmap_dbl(list(data, lead(data), CommonID),
                                          ~ {if (all(is.na(..2))) {
                                              return(NA)
                                          } else {                                    
                                              (filter(..2, ID %in% ..3) %$% mean(Pheno_scale)) -
                                                  (filter(..1, ID %in% ..3) %$% mean(Pheno_scale))
                                          }}),
               # And its standard error!
               DeltaCommonZbar_SE = pmap_dbl(list(data, lead(data), CommonID),
                                             ~ {if (all(is.na(..2))) {
                                                 return(NA)
                                             } else {                                    
                                                 sd(..2[match(..3, ..2[["ID"]]), ][["Pheno_scale"]] -
                                                    ..1[match(..3, ..1[["ID"]]), ][["Pheno_scale"]]) /
                                                   length(..3)
                                             }})) %>%
        # Now, we need to complete the missing years (to match thetas)
        complete(Year = full_seq(Year, 1))
    
    # Cleaning up the output
    deltas <-
        deltas %>%
        select(-data, -ID, -CommonID) %>%
        slice(-n())                         # Removing last line because it's useless
    
    return(deltas)
}

## Function to compute the row-wise correlations between two matrices
# Args: - mat1: a matrix
#       - mat2: a matrix of dim equal to mat1
# Value: a vector of correlation values
rowwise_correlate <- function(mat1, mat2, method = "pearson") {
    cor <- numeric(nrow(mat1))
    for (i in 1:nrow(mat1)) {
        cor[i] <- cor(mat1[i, ], mat2[i, ], method = method)
    }
    return(cor)
}


## -------------------------------- Mathematical functions

## Gaussian optimum model
# Args: - x: Phenological trait
#       - theta: Location of the optimum
#       - wmax: Maximum of fitness (at the optimum)
#       - omega: Width parameter of the optimum
# Value: Expected fitness at x with the passed arguments (type double, possibly vectorised)
opt_model <- function(x, theta, wmax, omega) {
    wmax * exp(-((x - theta)^2 / (2 * omega^2)))
}

## Gaussian optimum model (logit)
# Args: - x: Phenological trait
#       - theta: Location of the optimum
#       - wmax: Maximum of fitness (at the optimum)
#       - omega: Width parameter of the optimum
# Value: Expected fitness at x with the passed arguments (type double, possibly vectorised)
opt_model_logit <- function(x, theta, wmax, omega) {
    plogis(log(wmax) - ((x - theta)^2 / (2 * omega^2)))
}

## Exponential model
# Args: - x: Phenological trait
#       - a: Intercept of the model
#       - b: Slope of the model
# Value: Expected fitness at x with the passed arguments (type double, possibly vectorised)
exp_model <- function(x, a, b) {
    exp(a + b * z)
}

## Exponential model (logit)
# Args: - x: Phenological trait
#       - a: Intercept of the model
#       - b: Slope of the model
# Value: Expected fitness at x with the passed arguments (type double, possibly vectorised)
exp_model_logit <- function(x, a, b) {
    plogis(a + b * z)
}

## Efficient log(1 - exp(-x))
# Args: -x: a numeric vector or scalar
# Value: the accurate value(s) of log(1 - exp(-x))
log1mexp <- function(x) {
    bool <- x < log(2)
    x[bool] <- log(-expm1(-x[bool]))
    x[!bool] <- log1p(-exp(-x[!bool]))
}

## Density of zero-inflated poisson
# Args: - x: point to estimate the density at
#       - lambda: Rate parameter of the Poisson
#       - zi: Zero-inflated probability of the Poisson
# Value: density of the ZI-Poisson
dzipois <- function(x, lambda, zi, log = FALSE) {
    # This way of computing the ZI-Poisson MPF should be fairly optimised
    lik <- numeric(length(x))
    bool <- x == 0 
    # Ensure the dimensions of lambda and zi
    if (length(zi) != 1) {
        zi <- zi[bool]
    }
    if (length(lambda) == 1) {
        lambda <- lambda * (numeric(length(x)) + 1)
    }    
    lik[bool] <-
        matrixStats::rowLogSumExps(
            cbind(dbinom(1, 1, zi, log = TRUE),
                  dbinom(0, 1, zi, log = TRUE) + dpois(0, lambda[bool], log = TRUE)),
        )
    lik[!bool] <- dpois(x[!bool], lambda[!bool], log = TRUE)
    if (!log) { lik <- exp(lik) }
    return(lik)
}

## Generating zero-inflated poisson data
# Args: - N: Number of samples to generate
#       - lambda: Rate parameter of the Poisson
#       - zi: Zero-inflated probability of the Poisson
# Value: Random draws from a ZI-Poisson (type int of size N)
rzipois <- function(N, lambda, zi) {
    # Check lengths
    if (!(length(lambda) %in% c(1, N))) {
        stop("lambda should of length 1 or N!")
    }
    if (!(length(zi) %in% c(1, N))) {
        stop("zi should of length 1 or N!")
    }
    
    trial <- runif(N) < zi
    out <- rpois(N, lambda = lambda)
    out <- if_else(trial,
                   0L,
                   out)
    return(out)
}

## Density of zero-truncated poisson
# Args: - x: point to estimate the density at
#       - lambda: Rate parameter of the Poisson
# Value: density of the zero-truncated Poisson
dztpois <- function(x, lambda, log = FALSE) {
    lik <- dpois(x, lambda, log = TRUE) - log1mexp(lambda)
    if (!log) { lik <- exp(lik) }
    return(lik)
}


## Sampling from a truncated Poisson
# Args: - N: Number of samples to generate
#       - lambda: Rate parameter of the Poisson
# Value: Random draws from a zero-truncated-Poisson (type int of size N)
rztpois <- function(N, lambda) {
    qpois(runif(N, dpois(0, lambda), 1), lambda)
}

## Function to transform GLMM parameters into Gaussian optimum ones
# Args: -b0, b1, b2: Parameters (intercept, slope of 1st, slope of 2nd order) of a GLMM
# Value: a tibble containing Theta, Wmax and Omega
glmm_to_explicit <- function(b0, b1, b2, dist = "Poisson") {
    if (dist %in% c("poisson", "ZIPoisson", "zeroinflatedpoisson1")) {
        dist <- "Poisson"
    }
    if (dist %in% c("Binom", "Binomial", "binomial")) {
        dist <- "Binom"
    }
    
    if (dist == "Poisson") {
        tibble(Omega = sqrt(-1 / (2 * b2)),
               Theta = - b1 / (2 * b2),
               Wmax  = exp(b0 - b1^2/(4 * b2)))
    } else if (dist == "Binom") {
        stop("Not implemented.")
    }

}

## Function to transform Gaussian into GLMM parameters optimum ones
# Args: -omega, theta, wmax: Parameters of the Gaussian optimum
# Value: a tibble containing b0, b1, b2 (intercept, slope of 1st, slope of 2nd order) of a GLMM
explicit_to_glmm <- function(omega, theta, wmax) {
    tibble(B0   = log(wmax) - (theta^2 / (2 * omega^2)),
           B1   = theta / (omega^2),
           B2   = - 1 / (2 * omega^2))
}

## Function to transform Wmax estimates back to log-scale
# Args:  - w_mu: intercept of Wmax on its own scale
#        - w_sigma: standard-deviation for Wmax on its own scale
# Value: a tbl containing the values for the intercept and 
#        standard deviation of Wmax on the log scale
transform_Wlog <- function(w_mu, w_sigma) {
    log_w_mu <- log((w_mu^3) / (w_sigma^2 + w_mu^2))
    log_w_sigma <- 2 * log(((w_sigma^2)/(w_mu^2)) + 1)
    
    return(tibble(log_w_intercept = log_w_mu, log_w_sigma = log_w_sigma))
}

## Compute a p-value from a vector
# Args:  - vec: a sample from a distribution (numeric)
# Value: a p-value testing whether the vector is different from 0 or not
compute_pval <- function(vec) {
    2 * sum((sign(median(vec)) * vec) < 0) / length(vec)
}

## Compute the gradient for a logit GLM
# Args: - est: estimates of the logit GLM (intercept and slope), a numeric vector of length 2
#       - z: the phenotypic trait used to fit the model
# Value: returns the gradient beta corresponding to the logit model
compute_beta_logit <- function(est, z) {
    if (length(est) != 2) {
        stop("There should be exactly two estimates")
    }
    
    # Computing the fitness function
    fit  <- plogis(est[1] + est[2] * z)
    
    # Computing the gradient
    beta <- est[2] * (1 - (mean(fit^2) / mean(fit)))
    
    return(beta)
}

## --------------------------------- Running and wrangling the models

## Defining a safe version of sampling
safe_sampling <- safely(rstan::sampling, otherwise = NULL)


## Fitting HMC
# Args: - data: a tbl containing the necessary data
#       - model: a stan_model object
#       - init: a list of initial values
#       - id: a unique identifier of the current model
#       - nchains: number of HMC chains to run
#       - progress: will a progress bar be used (must be named pb and use "progress" pkg)
#       - control: list of control parameters to pass to "rstan::sampling"
#       - infile: should output be saved in a file (TRUE) or returned (FALSE)
# Value: A STAN model
fit_hmc <- function(data,
                    model,
                    id,
                    nchains     = 4,
                    progress    = FALSE,
                    control     = list(adapt_delta = 0.80, max_treedepth = 12),
                    infile      = TRUE,
                    noind       = FALSE) {
    
    if (progress) { pb$tick() }
    
    # Type of model
    type_year_model <-
        model@model_name %>%
        stringr::str_match("type([A-Z1]+(_withtrend)?)") %>%
        pluck(2)
    
    # Grouping factors
    J_year  <- data[["Year_fac"]] %>% as.integer()
    J_ind   <- data[["ID"]] %>% forcats::as_factor() %>% as.integer()
    
    # Create the inits for each chain
    # Need to account for ZIPoisson case
    if (stringr::str_detect(model@model_name, "ZIPoisson")) {
        keep_fit <- data[["Fitness"]] != 0
    } else {
        keep_fit <- rep(TRUE, nrow(data))
    }
    init_log_Wmax   <- with(data, log(mean(Fitness[keep_fit])))
    init_pzi        <- with(data, mean(Fitness == 0))
    Nyear           <- max(J_year)
    Nind            <- max(J_ind)
    init <- list(omega      = 1,
                 theta      = 0, 
                 log_Wmax   = init_log_Wmax,
                 intercept  = init_log_Wmax,
                 p_zi       = init_pzi,
                 phi        = 0,
                 trend      = 0,
                 slope      = 0,
                 t_intercept= 0,
                 t_year     = rep(0, Nyear),
                 sl_year    = rep(0, Nyear),
                 t_sigma    = 1,
                 sl_sigma   = 1,
                 log_w_intercept= init_log_Wmax,
                 log_w_year     = rep(init_log_Wmax, Nyear),
                 int_year       = rep(init_log_Wmax, Nyear),
                 log_w_sigma    = 0.5,
                 int_sigma      = 0.5,
                 ind_sigma      = 0.5,
                 log_w_ind      = rep(0, Nind),
                 int_ind        = rep(0, Nind),
                 log_w_extend   = rep(0, Nyear),
                 t_extend       = rep(0, Nyear),
                 int_extend     = rep(0, Nyear),
                 sl_extend      = rep(0, Nyear),
                 ind_extend     = rep(0, Nind))
    init_chains <- map(1:nchains, ~ init)
    
    # Formatting data for STAN
    if (type_year_model == "NULL") {
        data_stan <-
            list(Nobs   = nrow(data),
                 W      = data[["Fitness"]],
                 z      = data[["Pheno_scale"]],
                 Nind   = Nind,
                 J_ind  = J_ind)
    } else if (type_year_model == "FLAT") {
        data_stan <-
            list(Nobs   = nrow(data),
                 W      = data[["Fitness"]],
                 Nyear  = Nyear,
                 J_year = J_year,
                 Nind   = Nind,
                 J_ind  = J_ind)
    } else {
        data_stan <-
            list(Nobs   = nrow(data),
                 W      = data[["Fitness"]],
                 z      = data[["Pheno_scale"]],
                 Nyear  = Nyear,
                 J_year = J_year,
                 Nind   = Nind,
                 J_ind  = J_ind)
    }
    
    # Parameters to keep
    if (type_year_model == "FLAT") {
        params <- c("lp__",
                    "log_w_ind", "ind_sigma",
                    "log_w_intercept", "log_w_sigma", "log_w_year")
    } else if (type_year_model == "NULL") {
        params <- c("lp__", "omega", "theta", "log_Wmax",
                    "log_w_ind", "ind_sigma")
    } else if (type_year_model == "VARINT") {
        params <- c("lp__", "omega", "theta",
                    "log_w_ind", "ind_sigma",
                    "log_w_intercept", "log_w_sigma", "log_w_year")
    } else if (type_year_model == "IID") {
        params <- c("lp__", "omega", "c",
                    "log_w_ind", "ind_sigma",
                    "t_intercept", "t_sigma", "t_year",
                    "log_w_intercept", "log_w_sigma", "log_w_year")
    } else if (type_year_model == "AR1") {
        params <- c("lp__", "omega", "c", "phi",
                    "log_w_ind", "ind_sigma",
                    "t_intercept", "t_sigma", "t_year",
                    "log_w_intercept", "log_w_sigma", "log_w_year")
    } else if (type_year_model == "AR1_withtrend") {
        params <- c("lp__", "omega", "c", "phi", "trend",
                    "log_w_ind", "ind_sigma",
                    "t_intercept", "t_sigma", "t_year",
                    "log_w_intercept", "log_w_sigma", "log_w_year")
    } else if (type_year_model == "EXPNULL") {
        params <- c("lp__", "intercept", "slope",
                    "int_ind", "ind_sigma")
    } else if (type_year_model == "EXPVARINT") {
        params <- c("lp__",
                    "int_ind", "ind_sigma",
                    "intercept", "int_sigma", "int_year",
                    "slope")
    } else if (type_year_model == "EXPIID") {
        params <- c("lp__",
                    "int_ind", "ind_sigma",
                    "intercept", "int_sigma", "int_year",
                    "slope", "sl_sigma", "sl_year")
    } else if (type_year_model == "EXPAR1") {
        params <- c("lp__", "phi",
                    "int_ind", "ind_sigma",
                    "intercept", "int_sigma", "int_year",
                    "slope", "sl_sigma", "sl_year")
    } else {
        stop("Invalid type of year model!")
    }
      
    # Adding p_zi of using a ZIPoisson
    if (stringr::str_detect(model@model_name, "ZIPoisson")) {
        params <- c(params, "p_zi")
    }
    
    if (noind) {
        params <- params[!str_detect(params, "ind")]
    }
    
    if (progress) {
        refresh <- -1
        show_messages <- FALSE
    } else {
        refresh <- 100
        show_messages <- TRUE
    }
    
    # Running the model
    fit_with_errors <- 
        safe_sampling(
            object = model,
            data   = data_stan,
            init   = init_chains,
            warmup = 1000,
            thin   = 5,
            iter   = 3000,
            control= control,
            pars   = params,
            chain  = nchains,
            refresh        = refresh,
            show_messages  = show_messages,
            #diagnostic_file= diag_file,
            open_progress  = FALSE
        )
    
    if (infile) {
        saveRDS(fit_with_errors,
                file = stringr::str_glue("../Output/{id}_model.rds"),
                version = 2)
        
        return(NULL)
    } else {
        return(fit_with_errors)
    }
}

## Diagnose a STAN object
# Args: - model: a stanfit object
#       - parameters: a list of parameters to study
# Value: a list of diagnostic values
diagnose_stan <- function(model, parameters) {   
    out <-
        list(
            Divergence  = rstan::get_divergent_iterations(model) %>% sum(),
            Energy      = rstan::get_bfmi(model) %>% {. < 0.2} %>% sum()
        )
}

## Diagnose a MCMC.list object
# Args: - mcmc: a mcmc.list object (possibly works with a plain mcmc as well)
#       - parameters: a list of parameters to study
# Value: a list of diagnostic values
diagnose_mcmc <- function(model, parameters) {
    tmp <- summary(model, pars = parameters, use_cache = FALSE)[["summary"]]
    
    n_itt <- model@sim[["chains"]] *
             (model@sim[["iter"]] - model@sim[["warmup"]]) / model@sim[["thin"]]
    
    out <-
        list(
            Rhat        = tmp[ , "Rhat"],
            N_eff       = tmp[ , "n_eff"],
            N_eff_ratio = tmp[ , "n_eff"] / n_itt
        )
}

## Function to read models and compute all diagnostics
# Args: - id: a fit ID (chr)
#       - model_year: the year model (chr)
#       - dist: the distribution (chr)
# Value: a tbl containg the mcmc and the diagnostic tests
load_and_diagnose <- function(id, model_year, dist, progress = TRUE, trace = TRUE) {
    # Update the progress bar
    if (progress) {
      pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Transform into MCMC
    params <- get_params_to_keep(model_year, dist)
    mcmc   <- rstan::As.mcmc.list(model, pars = params)
    
    # Plot the traceplot
    if (trace) {
        suppressMessages(
            ggsave(plot       = rstan::traceplot(model, pars = params),
                   filename   = str_glue("../Diagnostic/TracePlots/{id}_trace.pdf"),
                   device     = cairo_pdf())
        )
        dev.off()
    }
    
    # Diagnoses and output
    out <- tibble(MCMC = list(mcmc),
                  Diag_STAN = list(diagnose_stan(model, params)),
                  Diag_MCMC = list(diagnose_mcmc(model, params)))
    
    # Garbage collection to be nice on the memory
    rm(model)
    rm(mcmc)
#     gc()
    
    return(out)
}

## Generate the posterior predictive check
# Args: - model: a stanfit object
#       - data: a tbl containing the data from which the fit was done
#       - id: an ID for naming the ECDF plot
# Value: a tbl of ppp-value (the function also generates a graphic if plot is TRUE)
generate_pred_check <- function(model, data, id = NULL, plot = TRUE) {
    
    if (is.null(id) & plot) {stop("'id' must be provided if plot is TRUE.")}
    
    # Getting and formatting the MCMC
    post <- extract(model)
    colnames(post[["log_w_ind"]])   <- unique(data[["ID"]])
    post[["log_w_ind"]]             <- post[["log_w_ind"]][ , data[["ID"]]]
    colnames(post[["t_year"]])      <- unique(data[["Year"]]) %>% sort()
    post[["t_year"]]                <- post[["t_year"]][ , as.character(data[["Year"]])]
    colnames(post[["w_year"]])      <- unique(data[["Year"]]) %>% sort()
    post[["w_year"]]                <- post[["w_year"]][ , as.character(data[["Year"]])]
    post[["omega"]]                 <- matrix(post[["omega"]], ncol = 1)[ , rep(1, nrow(data))]
    
    # Computing the latent variable
    latent <- post[["log_w_ind"]] - 
              ((data[["Pheno_scale"]] - post[["t_year"]])^2) /
              (2 *  post[["omega"]]^2)
    latent <- post[["w_year"]] * exp(latent)
    
    # Generating new data (note that this "transposes" iterations and data points)
    # If no zero in the data, then use a truncated Poisson
    if (min(data[["Fitness"]]) > 0) {
        yrep <- apply(latent, 1, function(lambda) {rztpois(length(lambda), lambda)})
    } else {
        yrep <- apply(latent, 1, function(lambda) {rpois(length(lambda), lambda)})
    }
    
    # Plotting a hundred ECDFs if asked for
    if (plot) {
        
        max <- min(max(data[["Fitness"]]), 5)
        
        # Formatting the simulated values
        tbl <- imap_dfc(0:max,
                        ~ colSums(yrep == .))
        colnames(tbl) <- str_c("C", 0:max)
        tbl <- mutate(tbl, Type = "Simulation")
        
        # Adding the data
        tbl <-
            bind_rows(
                tbl,
                add_column(
                    imap_dfc(0:max,
                         ~ sum(data[["Fitness"]] == .)) %>%
                        set_colnames(str_c("C", 0:max)),
                    Type = "Data"
                )
            )
        
        # Plotting everything
        tbl <- gather(tbl, "At", "Count", -Type)
        p <- ggplot() +
             geom_violin(data = tbl %>% filter(Type == "Simulation"),
                         mapping = aes(x = At, y = Count),
                         fill = "grey80") +
             geom_point(data = tbl %>% filter(Type == "Data"),
                        mapping = aes(x = At, y = Count),
                        colour = "red", size = 5)
        suppressMessages(
            ggsave(filename = str_glue("../Diagnostic/PPC/{id}_post_check.pdf"),
                   plot     = p,
                   device   = cairo_pdf())
        )
    }
    
    # Computing the ppp-values
    ppp_zero <- mean(colSums(yrep == 0) < sum(data[["Fitness"]] == 0))
    ppp_one  <- mean(colSums(yrep == 1) < sum(data[["Fitness"]] == 1))
    
    return(tibble(ppp_zero = ppp_zero, ppp_one = ppp_one))
}

## Function to read models and perform posterior predictive check
# Args: - id: a fit ID (chr)
#       - data: the corresponding dataset (a tbl)
# Value: a tbl containg the ppp-values for the model
load_and_ppc <- function(id, data, progress = TRUE, plot = TRUE) {
    # Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Performing the posterior predictive check
    ppp <- generate_pred_check(model, data, id, plot = plot)
    
    rm(model)
    
    return(ppp)
}


## Merging chains from two models
# Args: - mod1, mod2: two stanfit objects (two instances of the exact same run!)
# Value: a stanfit object with merged chains from the two objects
merge_models <- function(mod1, mod2) {
    
    # Getting the number of chains for each
    nchains1 <- mod1@sim[["chains"]]
    nchains2 <- mod2@sim[["chains"]]
    
    # Using mod1 as the "template"
    mod <- mod1
    for (i in 1:nchains2) {
        attr(mod2@sim[["samples"]][[i]], "args")[["chain_id"]] <- nchains1 + i
    }
    
    mod@sim[["samples"]] <- c(mod1@sim[["samples"]], mod2@sim[["samples"]])
    mod@sim[["chains"]] <- nchains1 + nchains2
    mod@sim[["n_save"]] <- c(mod1@sim[["n_save"]], mod2@sim[["n_save"]])
    mod@sim[["warmup2"]] <- c(mod1@sim[["warmup2"]], mod2@sim[["warmup2"]])
    mod@sim[["permutation"]] <- c(mod1@sim[["permutation"]], mod2@sim[["permutation"]])
    mod@inits <- c(mod1@inits, mod2@inits)
    mod@stan_args <- c(mod1@stan_args, mod2@stan_args)
    
    return(mod)
}

## Replacing chains of a model by chains from another (when only a few failed...)
# Args: - mod1, mod2: two stanfit objects (two instances of the exact same run!)
#       - chains1: chains to be replaced in mod1 (integer vector)
#       - chains2: chains to replace chains1 with from mod2 (integer vector, length of chains1)
# Value: a stanfit object corresponding to mod1 with replaced chains
replace_chains <- function(mod1, mod2, chains1, chains2) {
    
    # Some checks
    if (length(chains1) != length(chains2)) {
        stop("Variables chains1 and chains2 must be of same length")
    }
        
    # Replacing the chains
    mod1@sim$samples[chains1]       <- mod2@sim$samples[chains2]
    mod1@sim$permutation[chains1]   <- mod2@sim$permutation[chains2]
    
    # Fixing chains id
    for (i in 1:length(mod1@sim$samples)) {
        attr(mod1@sim$samples[[i]], "args")[["chain_id"]] <- i
    }
    
    return(mod1)
}

## Compute likelihood array
# Args: - mcmc: an mcmc array, obtained using rstan::extract
#       - model_name: "model_name" slot from the stanfit object
#       - data: the corresponding dataset used to fit the stanfit object
# Comment: Not using the stanfit object directly to be gentle on memory usage...
# Value: a likelihood array formatted for the `loo` package
compute_likarray <- function(mcmc, model_name, data) {
    
    ## Getting information about the model
    dist    <- str_split(model_name, "_") %>%
               chuck(1, 2) %>%
               str_remove("dist")
    modyear <- str_split(model_name, "_") %>%
               chuck(1, 3) %>%
               str_remove("type")
    niter   <- length(mcmc[["lp__"]])
    
    ## Setting up the parameters in a tbl
    # Grouping factors
    J_year  <- data[["Year_fac"]] %>% as.integer()
    J_ind   <- data[["ID"]] %>% forcats::as_factor() %>% as.integer()
    
    # tbl will be used to store the parameters for logLik
    tbl <- tibble(y = map(1:niter, ~ data[["Fitness"]]),
                  z = map(1:niter, ~ data[["Pheno_scale"]]))
    
    # If FLAT model, then phenology (z) is not a parameter for logLik, let's throw it away
    if (modyear == "FLAT") {
        tbl <- tbl %>% select(-z)
    }
    
    # Formatting the MCMC parameters into something usable
    omega   <- mcmc[["omega"]]          # NULL if missing
    zi      <- mcmc[["p_zi"]]
    # Theta
    if (modyear %in% c("IID", "AR1")) {
        theta <- mcmc[["t_year"]][ , J_year]
        dimnames(theta) <- NULL        # Removing the names to avoid some mess after
        # Place each iteration in the elements of a list
        theta <- map(1:nrow(theta), ~ theta[., ])
    } else {
        theta <- mcmc[["theta"]]    # NULL if missing
        names(theta) <- NULL        # Removing the names to avoid some mess after
    }
    # Wmax
    if (modyear %in% c("FLAT", "IID", "AR1", "VARINT")) {
        log_wmax <- mcmc[["log_w_year"]][ , J_year] + mcmc[["log_w_ind"]][ , J_ind]
        dimnames(log_wmax) <- NULL        # Removing the names to avoid some mess after
        # Place each iteration in the elements of a list
        log_wmax <- map(1:nrow(log_wmax), ~ log_wmax[., ])
    } else if (modyear == "NULL") {
        log_wmax <- matrix(rep(mcmc[["log_Wmax"]], nrow(data)),
                           nrow = nrow(mcmc[["log_w_ind"]])) +
                    mcmc[["log_w_ind"]][ , J_ind]
        dimnames(log_wmax) <- NULL         # Removing the names to avoid some mess after
        # Place each iteration in the elements of a list
        log_wmax <- map(1:nrow(log_wmax), ~ log_wmax[., ])
    } else {
        log_wmax <- NULL
    }
    # Intercept/Slope for EXP models
    if (modyear %in% c("EXPIID", "EXPAR1")) {
        intercept <- mcmc[["int_year"]][ , J_year] + mcmc[["int_ind"]][ , J_ind]
        dimnames(intercept) <- NULL        # Removing the names to avoid some mess after
        # Place each iteration in the elements of a list
        intercept <- map(1:nrow(intercept), ~ intercept[., ])
        
        slope <- mcmc[["sl_year"]][ , J_year]
        dimnames(slope) <- NULL        # Removing the names to avoid some mess after
        # Place each iteration in the elements of a list
        slope <- map(1:nrow(slope), ~ slope[., ])
    } else if (modyear == "EXPVARINT") {
        intercept <- mcmc[["int_year"]][ , J_year] + mcmc[["int_ind"]][ , J_ind]
        dimnames(intercept) <- NULL        # Removing the names to avoid some mess after
        # Place each iteration in the elements of a list
        intercept <- map(1:nrow(intercept), ~ intercept[., ])
        
        slope <- mcmc[["slope"]]
        names(slope) <- NULL        # Removing the names to avoid some mess after
    } else if (modyear == "EXPNULL") {
        intercept <- matrix(rep(mcmc[["intercept"]], nrow(data)),
                            nrow = nrow(mcmc[["int_ind"]])) + 
                     mcmc[["int_ind"]][ , J_ind]
        dimnames(intercept) <- NULL        # Removing the names to avoid some mess after
        # Place each iteration in the elements of a list
        intercept <- map(1:nrow(intercept), ~ intercept[., ])
        
        slope <- mcmc[["slope"]]
        names(slope) <- NULL        # Removing the names to avoid some mess after
    } else {
        intercept   <- mcmc[["intercept"]]      # NULL if missing
        slope       <- mcmc[["slope"]]          # NULL if missing
    }
    
    # Put everything in the tbl (NULL values will result in no added column)
    tbl <-
        tbl %>%
        mutate(omega    = omega,
               theta    = theta,
               log_wmax = log_wmax,
               intercept= intercept,
               slope    = slope,
               zi       = zi)
    
    # Constructing the likelihood function
    if (dist == "Binom" & modyear == "FLAT") {
        func_logLik <- function(y, log_wmax) {
            dbinom(y, 1, plogis(log_wmax), log = TRUE)
        }
    } else if (dist == "Binom" & (modyear %in% c("NULL", "VARINT", "IID", "AR1"))) {
        func_logLik <- function(y, z, omega, theta, log_wmax) {
            p <- plogis(log_wmax) * exp(-((z - theta)^2) / (2 * omega^2))
            dbinom(y, 1, p, log = TRUE)
        }
    } else if (dist == "Binom" & str_detect(modyear, "EXP")) {
        func_logLik <- function(y, z, intercept, slope) {
            p <- intercept + slope * z
            dbinom(y, 1, plogis(p), log = TRUE)
        }
    } else if (dist == "Poisson" & modyear == "FLAT") {
        func_logLik <- function(y, log_wmax) {
            dpois(y, exp(log_wmax), log = TRUE)
        }
    } else if (dist == "Poisson" & (modyear %in% c("NULL", "VARINT", "IID", "AR1"))) {
        func_logLik <- function(y, z, omega, theta, log_wmax) {
            l <- log_wmax - (((z - theta)^2) / (2 * omega^2))
            dpois(y, exp(l), log = TRUE)
        }
    } else if (dist == "Poisson" & str_detect(modyear, "EXP")) {
        func_logLik <- function(y, z, intercept, slope) {
            l <- intercept + slope * z
            dpois(y, exp(l), log = TRUE)
        }
    } else if (dist == "TruncPoisson" & modyear == "FLAT") {
        func_logLik <- function(y, log_wmax) {
            dztpois(y, exp(log_wmax), log = TRUE)
        }
    } else if (dist == "TruncPoisson" & (modyear %in% c("NULL", "VARINT", "IID", "AR1"))) {
        func_logLik <- function(y, z, omega, theta, log_wmax) {
            l <- log_wmax - (((z - theta)^2) / (2 * omega^2))
            dztpois(y, exp(l), log = TRUE)
        }
    } else if (dist == "TruncPoisson" & str_detect(modyear, "EXP")) {
        func_logLik <- function(y, z, intercept, slope) {
            l <- intercept + slope * z
            dztpois(y, exp(l), log = TRUE)
        }
    } else if (dist == "ZIPoisson" & modyear == "FLAT") {
        func_logLik <- function(y, log_wmax, zi) {
            dzipois(y, exp(log_wmax), zi, log = TRUE)
        }        
    } else if (dist == "ZIPoisson" & (modyear %in% c("NULL", "VARINT", "IID", "AR1"))) {
        func_logLik <- function(y, z, omega, theta, log_wmax, zi) {
            l <- log_wmax - (((z - theta)^2) / (2 * omega^2))
            dzipois(y, exp(l), zi, log = TRUE)
        }
    } else if (dist == "ZIPoisson" & str_detect(modyear, "EXP")) {
        func_logLik <- function(y, z, intercept, slope, zi) {
            l <- intercept + slope * z
            dzipois(y, exp(l), zi, log = TRUE)
        }
    } else {
        stop("No defined likelihood function!")
    }
    
    # Computing the logLik for each iteration, then format as a matrix
    out <- pmap(tbl, func_logLik) %>%
            {do.call("rbind", .)} %>%
            as.matrix()
    
    return(out)
}

## Function to read models and perform leave-one-out (LOO)
# Args: - id: a fit ID (chr)
#       - data: the corresponding dataset (a tbl)
# Value: a loo object
load_and_loo <- function(id, data, progress = TRUE) {
    # Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Compute log-lik array
    mcmc <- rstan::extract(model)
    model_name <- model@model_name
    loo <- compute_likarray(mcmc, model_name, data) %>%
           loo::loo()
    
    # Garbage collection to be nice on the memory
    rm(model)
    rm(mcmc)
    gc()
    
    return(loo)
}

## Compute predicted values
# Args: - mcmc: an mcmc array, obtained using rstan::extract
#       - model_name: "model_name" slot from the stanfit object
#       - data: the corresponding dataset used to fit the stanfit object
# Comment: Not using the stanfit object directly to be gentle on memory usage...
# Value: a vector of predicted values (same length as nrow(data))
compute_predict <- function(mcmc, model_name, data) {
    
    ## Getting information about the model
    dist    <- str_split(model_name, "_") %>%
               pluck(1, 2) %>%
               str_remove("dist")
    modyear <- str_split(model_name, "_") %>%
               pluck(1, 3) %>%
               str_remove("type")
    niter   <- length(mcmc[["lp__"]])
    
    ## Getting the posterior medians to work on
    #NOTE Removing the individual values
    mcmc <- mcmc[!str_detect(names(mcmc), "_ind$")]
    #FIXME Maybe better to predict values on mcmc chains directly, but much longer!
    medians <- map(mcmc, ~ apply(as.matrix(.), 2, median)) %>%
               as_tibble()
    
    # Formatting the medians into something usable
    medians <-
        medians %>%
        # Filtering stuff useless for prediction
        mutate(lp__             = NULL,
               p_zi             = NULL,
               phi              = NULL,
               phi_w            = NULL,
               phi_t            = NULL,
               rho              = NULL,
               c                = NULL,
               t_intercept      = NULL,
               t_sigma          = NULL,
               log_w_intercept  = NULL,
               log_w_sigma      = NULL,
               int_sigma        = NULL,
               sl_sigma         = NULL,
               ind_sigma        = NULL)
    if (colnames(medians) %>% str_detect("int_year") %>% any()) {
        medians <-
            medians %>%
            select(-intercept, -slope)
    }
    medians <-
        medians %>%
        rename_at(vars(contains("_year"), contains("Wmax")),
                  ~ recode(.,
                           t_year    = "theta",
                           log_w_year= "log_wmax",
                           int_year  = "intercept",
                           sl_year   = "slope",
                           log_Wmax  = "log_wmax"))  
    
    # If all parameters were constant, "medians" has now 1 row and needs to be fixed
    if (nrow(medians) == 1) {
        nyear <- data[["Year"]] %>% unique() %>% length()
        medians <-
            crossing(Year = 1:nyear, medians) %>%
            select(-Year)
    }
    
    ## Setting up the parameters for computation
    # Grouping factors
    J_year  <- data[["Year_fac"]] %>% as.integer()
    medians <- medians[J_year, ]
    
    ## Adding the covariate if needed
    if (modyear != "FLAT") {
        medians <-
            medians %>%
            add_column(z = data[["Pheno_scale"]], .before = 1)
    }
    
    ## Constructing the latent predictor function
    if (modyear == "FLAT") {
        func_latent <- function(log_wmax) {
            log_wmax
        }
    } else if ((modyear %in% c("NULL", "VARINT", "IID", "AR1"))) {
        if (dist == "Binom") {
            # Needs to directly define the model
            func_latent <- function(z, omega, theta, log_wmax) {
                plogis(log_wmax) * exp(- ((z - theta)^2) / (2 * omega^2))
            }
        } else {
            func_latent <- function(z, omega, theta, log_wmax) {
                log_wmax - ((z - theta)^2) / (2 * omega^2)
            }
        }
    } else if (str_detect(modyear, "EXP")) {
        func_latent <- function(z, intercept, slope) {
            intercept + slope * z
        }
    } else {
        stop("No defined latent function!")
    }
    
    ## Constructing the link function
    if (dist == "Binom") {
        if ((modyear %in% c("NULL", "VARINT", "IID", "AR1"))) {
            # Identity function as the "latent" function is complex
            func_link <- function(x) {x}
        } else {
            func_link <- plogis
        }
    } else {
        func_link <- exp
    }
    
    # Computing the logLik for each iteration, then format as a matrix
    out <- do.call("func_latent", medians) %>%
           func_link()
    
    return(out)
}

## Function to read models and compute predicted values
# Args: - id: a fit ID (chr)
#       - data: the corresponding dataset (a tbl)
# Value: a tbl containg the ppp-values for the model
load_and_predict <- function(id, data, progress = TRUE) {
    # Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Compute log-lik array
    mcmc <- rstan::extract(model)
    model_name <- model@model_name
    pred <- compute_predict(mcmc, model_name, data)
    
    return(pred)
}

## Compute point estimates of a model
compute_point_estimates <- function(mcmc) {
    
    # Need to be defensive in the case of missing values
    if (any(is.na(mcmc))) {
        out <- as.list(rep(NA, 5)) %>%
               set_names(c("Median", "Mode", "SE", "Low", "Up"))
    } else {
        out <- flatten(list(
            Median  = median(mcmc),
            Mode    = MCMCglmm::posterior.mode(coda::as.mcmc(mcmc)),
            SE      = sd(mcmc),
            as.list(coda::HPDinterval(coda::as.mcmc(mcmc))[1, ])
        )) %>% set_names(c("Median", "Mode", "SE", "Low", "Up"))
    }
    
    return(out)
}

## Generate a summary of the model estimates
# Args: - model: a stanfit object
# Value: a tbl of summary statistics for each posterior distribution
summarise_rstan <- function(model) {
    
    ## Getting information about the model
    model_name <- model@model_name
    dist    <- str_split(model_name, "_") %>%
               chuck(1, 2) %>%
               str_remove("dist")
    modyear <- str_split(model_name, "_") %>%
               chuck(1, 3) %>%
               str_remove("type")
    
    ## Extracting the parameters of interest
    params <- get_params_to_keep(modyear, dist)
    post <- rstan:::as.data.frame.stanfit(model, pars = params) %>%
            as_tibble()
    
    ## Computing the summarise statistics of interest
    out <-
        post %>%
        map_dfr(~ tibble(!!!compute_point_estimates(.)),
                .id = "Parameter")
    
    return(out)
}

## Generate a summary of the brms meta-analysis
# Args: - model: a brmsfit object from the "meta-analysis" step
#       - sqrt: was a sqrt transform used?
# Value: a tbl of point estimates and summary statistics for the parameters of the model
summarise_meta_brms <- function(model) {
    
    # Extracting the MCMC samples
    mcmc <- 
        as.data.frame(model) %>%
        as_tibble() %>%
        select(-lp__)
    
    # Computing the estimates
    if (model[["family"]][["link"]] == "softplus") {
        mcmc <-
            mcmc %>%
            mutate(b_TaxonMammal = b_Intercept + b_TaxonMammal,
                   v_tot         = select(., starts_with("sd"), "sigma") %>%
                                   mutate(across(everything(), raise_to_power, 2)) %>%
                                   rowSums(),
                   mu_Bird       = map2_dbl(b_Intercept, v_tot,
                                        ~  QGglmm::QGmean(mu     = .x,
                                                          var    = .y,
                                                          link.inv = function(x) {
                                                              log(exp(x) + 1)
                                                          })),
                   mu_Mammal     = map2_dbl(b_TaxonMammal, v_tot,
                                        ~  QGglmm::QGmean(mu     = .x,
                                                          var    = .y,
                                                          link.inv = function(x) {
                                                              log(exp(x) + 1)
                                                          })),)
    } else {
        mcmc <-
            mcmc %>%
            mutate(b_TaxonMammal = b_Intercept + b_TaxonMammal) %>%
            rename(mu_Bird   = "b_Intercept",
                   mu_Mammal = "b_TaxonMammal")
    }
    
    # Dropping everything else other than the meta-estimates
    mcmc <- select(mcmc, mu_Bird, mu_Mammal)
    
    # Summarising the mcmc
    out <-
        mcmc %>%
        map_dfr(~ tibble(!!!compute_point_estimates(.)),
                .id = "Parameter")
    
    return(out)
}

## Function to read models and summarise them
# Args: - id: a fit ID (chr)
# Value: a tbl containg the ppp-values for the model
load_and_summarise <- function(id, progress = TRUE) {
    
    # Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Summarising the parameters
    out <- summarise_rstan(model)
    
    return(out)
}

## Function to read models and get theta from them
# Args: - id: a fit ID (chr)
# Value: posterior distribution of theta BLUPs from the fit (an array)
load_and_get_theta <- function(id, progress = TRUE) {
    
    # Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Getting theta
    par <- if_else(str_detect(id, "NULL"), "theta", "t_year")
    theta <- rstan::extract(model, pars = par)[[par]]
    
    return(theta)
}

## Function to read models and get omega from them
# Args: - id: a fit ID (chr)
# Value: posterior distribution of theta BLUPs from the fit (an array)
load_and_get_omega <- function(id, progress = TRUE) {
    
    # Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Getting theta
    omega <- rstan::extract(model, pars = "omega")[["omega"]]
    
    return(omega)
}

## Function to read models and get omega from them
# Args: - id: a fit ID (chr)
# Value: posterior distribution of theta BLUPs from the fit (an array)
load_and_get_t_sigma <- function(id, progress = TRUE) {
    
    # Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Getting theta
    t_sigma <- rstan::extract(model, pars = "t_sigma")[["t_sigma"]]
    
    return(t_sigma)
}

## Function to read models and get slope (b) from them
# Args: - id: a fit ID (chr)
# Value: posterior distribution of theta BLUPs from the fit (an array)
load_and_get_estexp <- function(id, progress = TRUE) {
    
    # Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    # Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    
    # Getting theta
    if (str_detect(id, "EXPNULL")) {
        par <- c("intercept", "slope")
    } else {
        par <- c("int_year", "sl_year")
    }

    int     <- rstan::extract(model, pars = par[1])[[par[1]]]
    slope   <- rstan::extract(model, pars = par[2])[[par[2]]]
    
    return(list(Intercept = int, Slope = slope))
}

## Generate the fitness function from the model
# Args: - id: a fit ID (chr)
#       - dist: a distribution name (chr)
#       - mod_year: a latent model name (chr)
# Value: the fitness function from the model (a function)
load_and_generate_fitfunc <- function(id, dist, mod_year, progress = TRUE) {
    
    ## Update the progress bar
    if (progress) {
        pb$tick()  
    }
    
    ## Reading the fit
    model <- readRDS(str_glue("../Output/{id}_model.rds"))[["result"]]
    est   <- rstan::extract(model)
    
    ## Selecting the correct link function
    if (dist == "Binom") {
        link_func <- plogis
    } else {
        link_func <- exp
    }
    
    ## Generating the fitness function
    if (mod_year == "FLAT") {
        
        # Getting parameters
        log_w_year  <- matrixStats::colMedians(est[["log_w_year"]])
        
        # Defining fitness function
        fit_func <- function(x, year_index) {
            link_func(log_w_year[year_index])
        }
        
    } else if (mod_year == "NULL") {
        
        # Getting parameters
        log_Wmax    <- median(est[["log_Wmax"]])
        theta       <- median(est[["theta"]])
        omega       <- median(est[["omega"]])
        
        # Defining fitness function
        if (dist == "Binom") {
            # Needs to directly define the model
            fit_func <- function(x, year_index) {
                link_func(log_Wmax) * exp(-((x - theta)^2) / (2 * omega^2))
            }
        } else {
            fit_func <- function(x, year_index) {
                link_func(log_Wmax - ((x - theta)^2 / (2 * omega^2)))
            }
        }
        
    } else if (mod_year == "EXPNULL") {
        
        # Getting parameters
        intercept   <- median(est[["intercept"]])
        slope       <- median(est[["slope"]])
        
        # Defining fitness function
        fit_func <- function(x, year_index) {
            link_func(intercept + slope * x)
        }
        
    } else if (mod_year == "VARINT") {
        
        # Getting parameters
        log_w_year  <- matrixStats::colMedians(est[["log_w_year"]])
        theta       <- median(est[["theta"]])
        omega       <- median(est[["omega"]])
        
        # Defining fitness function
        if (dist == "Binom") {
            # Needs to directly define the model
            fit_func <- function(x, year_index) {
                link_func(log_w_year[year_index]) * exp(-((x - theta)^2) / (2 * omega^2))
            }
        } else {
            fit_func <- function(x, year_index) {
                link_func(log_w_year[year_index] - ((x - theta)^2 / (2 * omega^2)))
            }
        }
        
    } else if (mod_year == "EXPVARINT") {
        
        # Getting parameters
        int_year    <- matrixStats::colMedians(est[["int_year"]])
        slope       <- median(est[["slope"]])
        
        # Defining fitness function
        fit_func <- function(x, year_index) {
            link_func(int_year[year_index] + slope * x)
        }
    } else if (mod_year %in% c("IID", "AR1")) {
        
        # Getting parameters
        log_w_year  <- matrixStats::colMedians(est[["log_w_year"]])
        theta_year  <- matrixStats::colMedians(est[["t_year"]])
        omega       <- median(est[["omega"]])
        
        # Defining fitness function
        if (dist == "Binom") {
            # Needs to directly define the model
            fit_func <- function(x, year_index) {
                link_func(log_w_year[year_index]) * exp(-((x - theta_year[year_index])^2) / (2 * omega^2))
            }
        } else {
            fit_func <- function(x, year_index) {
                link_func(log_w_year[year_index] - ((x - theta_year[year_index])^2 / (2 * omega^2)))
            }
        }
        
    } else if (mod_year %in% c("EXPIID", "EXPAR1")) {
        
        # Getting parameters
        int_year    <- matrixStats::colMedians(est[["int_year"]])
        sl_year     <- matrixStats::colMedians(est[["sl_year"]])
        
        # Defining fitness function
        fit_func <- function(x, year_index) {
            link_func(int_year[year_index] + sl_year[year_index] * x)
        }
    }
    
    rm(model)
    rm(est)
    
    return(fit_func)
}

## --------------------------------- Graphical functions

## Function to plot the data and the fitted values
# Args: - tbl: a tbl containing both original data and fitted values
# Value: a ggplot graph displaying the data and fitted values
plot_fitted <- function(tbl, with_smooth = FALSE, se_smooth = TRUE, scales = "fixed") {
    p <- 
        ggplot(tbl) +
        geom_count(aes(x = Pheno, y = Fitness), alpha = 0.5) +
        geom_line(aes(x = Pheno, y = Fitted),
                  colour    = "#0055ff",
                  size      = 1) +
        facet_wrap(~ Year_fac, scales = scales) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        labs(x = "Phenological trait", y = "Fitness trait") +
        scale_size_area(name = "Counts")
    if (with_smooth) {
        p <- p +
             geom_smooth(aes(x = Pheno, y = Fitness),
                         colour = "#aa000088",
                         size   = 0.8,
                         se     = se_smooth)
    }
    return(p)
}

## Function to generate a dataframe to plot fitness curve
# Args: - tbl: a tbl containing the data (tbl)
#       - fitfunc: the fitness function we want to plot (function)
# Value: a "gridded" tbl for plotting of the curve
construct_curve_df <- function(tbl, fitfunc) {
    # Getting the current year for transforming Pheno to a numerical
    current_year <- lubridate::year(lubridate::today())
    
    # Create a curve dataframe to plot fitness function
    curve_df <-
        crossing(Pheno = seq(min(tbl[["Pheno"]]),
                             max(tbl[["Pheno"]]),
                             length.out = 100),
                 Year_index = 1:nlevels(tbl[["Year_fac"]])) %>%
        mutate(Pheno_num   = {Pheno - lubridate::ymd(stringr::str_glue("{current_year}-01-01"))} %>%
                             as.double(),
               Pheno_scale = (Pheno_num - mean(tbl[["Pheno_num"]], na.rm = TRUE)) /
                             unique(tbl[["SD_Within"]]),
               Y = map2_dbl(Pheno_scale, Year_index,
                            ~ fitfunc(.x, .y)),
               Year_fac = levels(tbl[["Year_fac"]])[Year_index] %>%
                          as.factor())
    
    return(curve_df)
}


## Function to plot the data and the fitness function values
# Args: - tbl: a tbl containing both original data and fitted values
# Value: a ggplot graph displaying the data and fitted values
plot_fitfunc <- function(curve_df, fitfunc, scales = "fixed") {
    
    # Now creating the plot
    p <- 
        ggplot(curve_df) +
        geom_line(mapping   = aes(x      = Pheno,
                                  y      = Y,
                                  colour = Year_index,
                                  group  = Year_fac),
                  size      = 1,
                  alpha     = 0.5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        labs(x = "Phenological trait", y = "Fitness trait") +
        scale_size_area(name = "Counts") +
        scale_y_continuous(limits = c(0, NA)) +
        scale_colour_gradient(low = "#0055ff", high = "#ffd500", guide = "none")
    return(p)
}

## Function to animate the data and the fitness function values
# Args: - tbl: a tbl containing both original data and fitted values
# Value: a ggplot graph displaying the data and fitted values
animate_fitfunc <- function(tbl, fitfunc, scales = "fixed") {
    # Constructing the curve df to plot the fitness function
    curve_df <- construct_curve_df(tbl, fitfunc)
    
    # Now creating the animation
    p_anim <- 
        ggplot(tbl) +
        geom_count(aes(x = Pheno, y = Fitness, group = Year_fac), alpha = 0.5) +
        geom_line(data      = curve_df,
                  mapping   = aes(x = Pheno, y = Y),
                  colour    = "#0055ff",
                  size      = 1) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        labs(x = "Phenological trait", y = "Fitness trait") +
        scale_size_area(name = "Counts") +
        gganimate::transition_states(Year_fac,
                                     state_length = 2,
                                     transition_length = 1) +
        labs(title = 'Year: {closest_state}') + 
        gganimate::enter_fade() + gganimate::enter_grow() +
        gganimate::exit_fade() + gganimate::exit_shrink() +
        gganimate::ease_aes('cubic-in-out')
    return(p_anim)
}
