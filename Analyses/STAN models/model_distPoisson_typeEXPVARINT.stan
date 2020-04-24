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
    real intercept;             // Latent grand intercept
    real slope;                 // Latent slope
    
    real<lower=0> int_sigma;    // Sigma of intercept random effect
    vector[Nyear] int_extend;   // Extended parameters for intercept
    
    real<lower=0> ind_sigma;    // Sigma of individuals effects
    vector[Nind] ind_extend;    // Individual effects
}
transformed parameters {
    vector[Nyear] int_year;     // Effect of years on intercept
    vector[Nind] int_ind;       // Effect of individual on log_Wmax
    
    // Extended parameters prior
    int_year = intercept + int_extend * int_sigma;
    
    // Extended parameters for individual effects
    int_ind = ind_extend * ind_sigma;
} 
model {
    
    // Latent variable
    vector[Nobs] latent;
    
    for (n in 1:Nobs) {
        latent[n] = int_ind[J_ind[n]] + int_year[J_year[n]] + slope * z[n];
    }
    
    // Likelihood
    W ~ poisson_log(latent);
    
    // Priors
    intercept ~ normal(0,1000);
    slope  ~ normal(0,1000);
    int_extend ~ normal(0,1);
    int_sigma ~ normal(0,1);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
}
generated quantities { 
}
