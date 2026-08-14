
// input data
data {
  int<lower=0> N;
  array[N] int oocysts;
  array[N] int spz_total;
}

// parameters accepted by the model
parameters {
  real<lower=0> beta;
  real<lower=0> alpha;
  real<lower=0> phi; // negative binomial dispersion par
}

// transformed parameters
transformed parameters {
  vector[N] mu = (to_vector(oocysts) * beta) ./ (1 + to_vector(oocysts) / alpha);
  // above uses vectorized form of the oocysts assay which is integers, allows for computations
}

// model to be estimated
model {
  spz_total ~ neg_binomial_2(mu, phi);
  beta ~ normal(4000, 1500);
  alpha ~ normal(25, 5);
  phi ~ exponential(0.01);
}

//
generated quantities{
    // simulations for PPC
    array[N] int sim_values_tot_spz = neg_binomial_2_rng(mu, phi);

    // log likelihood for comparison to other models
    vector[N] log_lik;
    for ( i in 1:N ) log_lik[i] = neg_binomial_2_lpmf( spz_total[i] | mu[i], phi );


}
