// STAN code
functions { 
} 
data {
    int<lower=1> Nobs;                          // Total number of observations
    int<lower=0> W[Nobs];                       // Fitness variable
    vector[Nobs] z;                             // Phenotype variable
    
    int<lower=1> Nind;                          // Total number of years
    int<lower=1,upper=Nind> J_ind[Nobs];        // Grouping indices for ID
}
parameters {
    real intercept;             // Latent intercept
    real slope;                 // Latent slope
    
    real<lower=0> ind_sigma;    // Sigma of individuals effects
    vector[Nind] ind_extend;    // Individual effects
}
transformed parameters {
    vector[Nind] int_ind;         // Effect of individual on log_Wmax
    
    // Extended parameters for individual effects
    int_ind = ind_extend * ind_sigma;
} 
model {
    
    // Latent variable
    vector[Nobs] latent = intercept + slope * z;
    
    for (n in 1:Nobs) {
        latent[n] += int_ind[J_ind[n]];
    }
    
    // Likelihood
    W ~ poisson_log(latent);
    for (n in 1:Nobs) {
	// From brms hurdle Poisson function
        target += log1m_exp(-exp(latent[n]));
    }
    
    // Priors
    intercept ~ normal(0,1000);
    slope  ~ normal(0,1000);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
}
generated quantities { 
}
