// Model of oocyst counts from control group children

// input data
data {
  int<lower=0> N;
  array[N] int oocyst_count;
}

// parameters accepted by the model
parameters {
  real<lower=0> mu; // mean of neg bin distribution
  real<lower=0> phi; // negative binomial dispersion par
}

// model to be estimated
model {
  oocyst_count ~ neg_binomial_2(mu, phi) T[1 , 99];
  mu ~ normal(10, 20); // approx because of mean and sd (oo_pos$oocyst_count)
  phi ~ exponential(0.01); // unsure, just basic exp
}

//
generated quantities{
    // simulations for PPC (removing the zeros and re-drawing)
    array[N] int sim_oocyst_counts;
    for ( i in 1:N ) sim_oocyst_counts[i] = neg_binomial_2_rng(mu, phi);

    // log likelihood for comparison to other models
    vector[N] log_lik;
    for ( i in 1:N ) log_lik[i] = neg_binomial_2_lpmf( oocyst_count[i] | mu, phi ) -
                                  log1m(neg_binomial_2_cdf( 0 | mu, phi ));
    // subtracting log(1 - P(X=0)) from the neg bin likelihood to normalize it based on zero-trncation (same as 1/1-dnbinom(0,mu, phi))
    // and subtracting log(1 - P(1 <))

}
// truncation for log-lik computes the log-density of the truncated distribution evaluated at the observed data points
