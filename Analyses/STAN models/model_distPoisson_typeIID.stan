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

    real t_intercept;                   // Intercept of theta
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
    vector[Nyear] t_year;           // Effect of years on theta
    vector[Nyear] log_w_year;       // Effect of years on log_Wmax
    vector[Nind] log_w_ind;         // Effect of individual on log_Wmax
    
    // Index theta_sigma to omega value
    t_sigma = c * omega;
    
    // Extended parameters prior
    log_w_year = log_w_intercept + log_w_extend * log_w_sigma;
    t_year = t_intercept + t_extend * t_sigma;
    
    // Extended parameters for individual effects
    log_w_ind = ind_extend * ind_sigma;
}  
model {
    
    // Latent variable
    vector[Nobs] latent;
    
    for (n in 1:Nobs) {
        latent[n] = log_w_ind[J_ind[n]] + log_w_year[J_year[n]] -
                    ((z[n] - t_year[J_year[n]]) * (z[n] - t_year[J_year[n]]) / (2 * omega * omega));
    }
    
    // Likelihood
    W ~ poisson_log(latent);
    
    // Priors
    omega ~ gamma(3.36,0.78);
    c ~ exponential(1);
    t_intercept ~ normal(0,1000);
    t_extend ~ normal(0,1);
    log_w_intercept ~ normal(0,1000);
    log_w_sigma ~ normal(0,1);
    log_w_extend ~ normal(0,1);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
}
generated quantities { 
}
