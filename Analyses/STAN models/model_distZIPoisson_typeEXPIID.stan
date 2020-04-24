// STAN code
functions { 
    // Zero-inflated Poisson from brms package
    /* zero-inflated poisson log-PDF of a single response 
     * Args: 
     *   y: the response value 
     *   lambda: mean parameter of the poisson distribution
     *   zi: zero-inflation probability
     * Returns:  
     *   a scalar to be added to the log posterior 
     */ 
    real zero_inflated_poisson_lpmf(int y, real lambda, real zi) { 
        if (y == 0) { 
            return log_sum_exp(bernoulli_lpmf(1 | zi), 
                               bernoulli_lpmf(0 | zi) + 
                               poisson_lpmf(0 | lambda)); 
        } else { 
            return bernoulli_lpmf(0 | zi) +  
            poisson_lpmf(y | lambda); 
        } 
    }
    /* zero-inflated poisson log-PDF of a single response 
     * logit parameterization of the zero-inflation part
     * Args: 
     *   y: the response value 
     *   lambda: mean parameter of the poisson distribution
     *   zi: linear predictor for zero-inflation part 
     * Returns:  
     *   a scalar to be added to the log posterior 
     */ 
    real zero_inflated_poisson_logit_lpmf(int y, real lambda, real zi) { 
        if (y == 0) { 
            return log_sum_exp(bernoulli_logit_lpmf(1 | zi), 
                               bernoulli_logit_lpmf(0 | zi) + 
                               poisson_lpmf(0 | lambda)); 
        } else { 
            return bernoulli_logit_lpmf(0 | zi) +  
            poisson_lpmf(y | lambda); 
        } 
    }
    /* zero-inflated poisson log-PDF of a single response
     * log parameterization for the poisson part
     * Args: 
     *   y: the response value 
     *   eta: linear predictor for poisson distribution
     *   zi: zero-inflation probability
     * Returns:  
     *   a scalar to be added to the log posterior 
     */ 
    real zero_inflated_poisson_log_lpmf(int y, real eta, real zi) { 
        if (y == 0) { 
            return log_sum_exp(bernoulli_lpmf(1 | zi), 
                               bernoulli_lpmf(0 | zi) + 
                               poisson_log_lpmf(0 | eta)); 
        } else { 
            return bernoulli_lpmf(0 | zi) +  
            poisson_log_lpmf(y | eta); 
        } 
    }
    /* zero-inflated poisson log-PDF of a single response 
     * log parameterization for the poisson part
     * logit parameterization of the zero-inflation part
     * Args: 
     *   y: the response value 
     *   eta: linear predictor for poisson distribution
     *   zi: linear predictor for zero-inflation part 
     * Returns:  
     *   a scalar to be added to the log posterior 
     */ 
    real zero_inflated_poisson_log_logit_lpmf(int y, real eta, real zi) { 
        if (y == 0) { 
            return log_sum_exp(bernoulli_logit_lpmf(1 | zi), 
                               bernoulli_logit_lpmf(0 | zi) + 
                               poisson_log_lpmf(0 | eta)); 
        } else { 
            return bernoulli_logit_lpmf(0 | zi) +  
            poisson_log_lpmf(y | eta); 
        } 
    }
    // zero-inflated poisson log-CCDF and log-CDF functions
    real zero_inflated_poisson_lccdf(int y, real lambda, real zi) { 
        return bernoulli_lpmf(0 | zi) + poisson_lccdf(y | lambda); 
    }
    real zero_inflated_poisson_lcdf(int y, real lambda, real zi) { 
        return log1m_exp(zero_inflated_poisson_lccdf(y | lambda, zi));
    }
} 
data {
    int<lower=1> Nobs;                          // Total number of observations
    int<lower=0> W[Nobs];                       // Fitness variable
    vector[Nobs] z;                             // Phenotype variable
    
    int<lower=1> Nyear;                         // Total number of years
    int<lower=1,upper=Nyear> J_year[Nobs];      // Grouping indices for years
    
    int<lower=1> Nind;                          // Total number of years
    int<lower=1,upper=Nind> J_ind[Nobs];        // Grouping indices for ID
} 
parameters {
    real intercept;                     // Latent grand intercept
    real slope;                         // Latent grand slope
    
    real<lower=0> int_sigma;            // Sigma of intercept random effect
    vector[Nyear] int_extend;           // Extended parameters for intercept
    
    real<lower=0> sl_sigma;             // Sigma of slope random effect
    vector[Nyear] sl_extend;            // Extended parameters for slope
    
    real<lower=0> ind_sigma;            // Sigma of individuals effects
    vector[Nind] ind_extend;            // Individual effects
    
    real<lower=0,upper=1> p_zi; // Probability of zero-inflation
}
transformed parameters { 
    vector[Nyear] int_year;     // Effect of years on intercept
    vector[Nyear] sl_year;      // Effect of years on slope
    vector[Nind] int_ind;       // Effect of individual on log_Wmax
    
    
    // Extended parameters prior
    int_year = intercept + int_extend * int_sigma;                  
    sl_year = slope + sl_extend * sl_sigma;
    
    // Extended parameters for individual effects
    int_ind = ind_extend * ind_sigma;
} 
model {
    
    // Latent variable
    vector[Nobs] latent;
    
    for (n in 1:Nobs) {
        latent[n] = int_ind[J_ind[n]] + int_year[J_year[n]] + sl_year[J_year[n]] * z[n];
    }
    
    // Likelihood
    for (n in 1:Nobs) {
        W[n] ~ zero_inflated_poisson_log_lpmf(latent[n], p_zi);
    }
    
    // Priors
    intercept ~ normal(0,1000);
    slope ~ normal(0,1000);
    int_extend ~ normal(0,1);
    int_sigma ~ normal(0,1);
    sl_extend ~ normal(0,1);
    sl_sigma ~ normal(0,1);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
    p_zi ~ uniform(0,1);
}
generated quantities { 
}
