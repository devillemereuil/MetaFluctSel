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
    real<lower=0> omega;            // Peak Width
    real theta;                     // Values of the optimum parameter
    
    real log_w_intercept;           // Intercept of log(Wmax)
    real<lower=0> log_w_sigma;      // Sigma of year effects on log(Wmax)
    vector[Nyear] log_w_extend;     // Extended parameters for log(Wmax)
    
    real<lower=0> ind_sigma;        // Sigma of individuals effects
    vector[Nind] ind_extend;        // Individual effects
    
    real<lower=0,upper=1> p_zi; // Probability of zero-inflation
}
transformed parameters {
    vector[Nyear] log_w_year;       // Effect of years on log_Wmax
    vector[Nind] log_w_ind;         // Effect of individual on log_Wmax
    
    // Extended parameters prior
    log_w_year = log_w_intercept + log_w_extend * log_w_sigma;
    
    // Extended parameters for individual effects
    log_w_ind = ind_extend * ind_sigma;
} 
model {
    
    // Latent variable
    vector[Nobs] latent;
    
    for (n in 1:Nobs) {
        latent[n] = log_w_ind[J_ind[n]] + log_w_year[J_year[n]] -
                    ((z[n] - theta) * (z[n] - theta) / (2 * omega * omega));
    }
    
    // Likelihood
    for (n in 1:Nobs) {
        W[n] ~ zero_inflated_poisson_log_lpmf(latent[n], p_zi);
    }
    
    // Priors
    omega ~ gamma(3.36,0.78);
    theta ~ normal(0,1000);
    log_w_intercept ~ normal(0,1000);
    log_w_sigma ~ normal(0,1);
    log_w_extend ~ normal(0,1);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
    p_zi ~ uniform(0,1);
}
generated quantities { 
}
