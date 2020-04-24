// STAN code
functions { 
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
    real<lower=0> omega;                // Peak width
    real<lower=-0.99,upper=0.99> phi;   // AR1 slope

    real process_t_mean;                // AR1 process mean for theta
    real<lower=0> c;                    // Proportionality coefficient to omega
    vector[Nyear] t_extend;             // Extended parameters for theta
    
    real log_w_intercept;               // Intercept of log(Wmax)
    real<lower=0> log_w_sigma;          // Sigma of year effects on log(Wmax)
    vector[Nyear] log_w_extend;         // Extended parameters for log(Wmax)
    
    real<lower=0> ind_sigma;            // Sigma of individuals effects
    vector[Nind] ind_extend;            // Individual effects
}
transformed parameters { 
    real<lower=0> t_sigma;          // Sigma of year effects on theta
    real<lower=0> process_t_sigma;  // AR1 process sd for theta
    real t_intercept;               // Intercept of theta
        
    vector[Nyear] t_year;           // Effect of years on theta
    vector[Nyear] log_w_year;       // Effect of years on log(Wmax)
    vector[Nind] log_w_ind;         // Effect of individual on log(Wmax)
    
    // Index process_theta_sigma to omega value
    process_t_sigma = c * omega;
    
    // Compute stationary mean and variance
    t_intercept = process_t_mean / (1 - phi);
    t_sigma = process_t_sigma / sqrt(1 - phi^2);
    
    // Auto-correlated process for theta
    t_year[1] = t_intercept + t_extend[1] * t_sigma;
    
    for (i in 2:Nyear) {
        t_year[i] = process_t_mean + 
                    phi * t_year[i-1] + 
                    t_extend[i] * process_t_sigma;
    }
    
    // Extended parameters for log(Wmax)
    log_w_year = log_w_intercept + log_w_extend * log_w_sigma;
    log_w_ind = ind_extend * ind_sigma;
}  
model {
    
    // Latent variable
    vector[Nobs] latent;
    
    for (n in 1:Nobs) {
        latent[n] = inv_logit(log_w_ind[J_ind[n]] + log_w_year[J_year[n]]) *
                    exp(-((z[n] - t_year[J_year[n]]) * (z[n] - t_year[J_year[n]]) / (2 * omega * omega)));
    }
    
    // Likelihood
    W ~ bernoulli(latent);
    
    // Priors
    omega ~ gamma(3.36,0.78);
    c ~ exponential(1);
    process_t_mean ~ normal(0,1000);
    t_extend ~ normal(0,1);
    log_w_intercept ~ normal(0,1);
    log_w_sigma ~ normal(0,1);
    log_w_extend ~ normal(0,1);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
    phi ~ uniform(-1,1); 
}
generated quantities { 
}
