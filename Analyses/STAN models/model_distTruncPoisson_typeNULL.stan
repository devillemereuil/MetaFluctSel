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
    real<lower=0> omega;            // Peak Width
    real theta;                     // Values of the optimum parameter
    real<upper=3> log_Wmax;         // Values of the max fitness parameter
    
    real<lower=0> ind_sigma;        // Sigma of individuals effects
    vector[Nind] ind_extend;        // Individual effects
}
transformed parameters {
    vector[Nind] log_w_ind;         // Effect of individual on log_Wmax
    
    // Extended parameters for individual effects
    log_w_ind = ind_extend * ind_sigma;
} 
model {
    
    // Latent variable
    vector[Nobs] latent = log_Wmax - ((z - theta) .* (z - theta)) ./ (2 * omega * omega);
    
    for (n in 1:Nobs) {
        latent[n] += log_w_ind[J_ind[n]];
    }
    
    // Likelihood
    W ~ poisson_log(latent);
    for (n in 1:Nobs) {
	// From brms hurdle Poisson function
        target += log1m_exp(-exp(latent[n]));
    }
    
    // Priors
    omega ~ gamma(3.36,0.78);
    theta ~ normal(0,1000);
    log_Wmax  ~ normal(0,1000);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
}
generated quantities { 
}
