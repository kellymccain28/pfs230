
// input data
data {
  int<lower=0> N;
  array[N] int oocysts;
  array[N] int spz_total;
}

// parameters accepted by the model
parameters {
  real<lower=0> mu; // r - constant value of total spz
  real<lower=0> phi; // negative binomial dispersion par
}

// model to be estimated
model {
  spz_total ~ neg_binomial_2(mu, phi);
  mu ~ lognormal(11, 1.1); // from mean(log(d1_spz_oocysts$spz_total)) and
  phi ~ exponential(0.01);
}

//
generated quantities{
    // simulations for PPC
    array[N] int sim_values_tot_spz;
    for ( i in 1:N ) sim_values_tot_spz[i] = neg_binomial_2_rng(mu, phi);

    // log likelihood for comparison to other models
    vector[N] log_lik;
    for ( i in 1:N ) log_lik[i] = neg_binomial_2_lpmf( spz_total[i] | mu, phi );


}
