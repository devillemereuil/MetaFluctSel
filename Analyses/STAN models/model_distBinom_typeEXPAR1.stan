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
    real<lower=-0.99,upper=0.99> phi;   // AR1 slope
    
    real intercept;                     // Grand intercept
    real<lower=0> int_sigma;            // Sigma of year effects on intercept
    vector[Nyear] int_extend;           // Extended parameters for intercept

    real process_sl_mean;               // AR1 process mean for the slope
    real<lower=0> process_sl_sigma;     // AR1 process sd for the slope
    vector[Nyear] sl_extend;            // Extended parameters for slope
    
    real<lower=0> ind_sigma;            // Sigma of individuals effects
    vector[Nind] ind_extend;            // Individual effects
}
transformed parameters { 
    real slope;                         // Latent grand slope
    real<lower=0> sl_sigma;             // Sigma of slope random effect
    
    vector[Nyear] int_year;         // Effect of years on intercept
    vector[Nyear] sl_year;          // Effect of years on slope
    vector[Nind] int_ind;           // Effect of individual on log_Wmax
        
    // Compute stationary mean and variance
    slope = process_sl_mean / (1 - phi);
    sl_sigma = process_sl_sigma / sqrt(1 - phi^2);
    
    // Auto- and cross-correlation process for intercept and slope
    sl_year[1] = slope + sl_extend[1] * sl_sigma;
    
    for (i in 2:Nyear) {
        sl_year[i] = process_sl_mean + 
                     phi * sl_year[i-1] + 
                     sl_extend[i] * process_sl_sigma;
    }
    
    // Extended parameters for the intercept
    int_year = intercept + int_extend * int_sigma;
    int_ind = ind_extend * ind_sigma;
} 
model {
    
    // Latent variable
    vector[Nobs] latent;
    
    for (n in 1:Nobs) {
        latent[n] = inv_logit(int_ind[J_ind[n]] + int_year[J_year[n]] + sl_year[J_year[n]] * z[n]);
    }
    
    // Likelihood
    W ~ bernoulli(latent);
    
    // Priors
    intercept ~ normal(0,1);
    int_sigma ~ normal(0,1);
    int_extend ~ normal(0,1);
    process_sl_mean ~ normal(0,1);
    sl_extend ~ normal(0,1);
    process_sl_sigma ~ normal(0,1);
    ind_extend ~ normal(0,1);
    ind_sigma ~ normal(0,1);
    phi ~ uniform(-1,1); 
}
generated quantities { 
}
