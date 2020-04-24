// STAN code
functions { 
} 
data {
    int<lower=1> Nobs;                          // Total number of observations
    int<lower=0> W[Nobs];                       // Fitness variable
    
    int<lower=1> Nyear;                         // Total number of years
    int<lower=1,upper=Nyear> J_year[Nobs];      // Grouping indices for years
    
    int<lower=1> Nind;                          // Total number of years
    int<lower=1,upper=Nind> J_ind[Nobs];        // Grouping indices for ID
} 
parameters {
    real log_w_intercept;           // Intercept of log_Wmax
    real<lower=0> log_w_sigma;      // Sigma of year effects on log_Wmax
    vector[Nyear] log_w_extend;     // Extended parameters for log(Wmax)

    real<lower=0> ind_sigma;        // Sigma of individuals effects
    vector[Nind] ind_extend;        // Individual effects
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
        latent[n] = log_w_ind[J_ind[n]] + log_w_year[J_year[n]];
    }
    
    // Likelihood
    W ~ poisson_log(latent);
    
    // Priors
    log_w_intercept ~ normal(0,1000);
    log_w_sigma ~ normal(0,1);
    log_w_extend ~ normal(0,1);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
}
generated quantities { 
}
