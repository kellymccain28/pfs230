// Model of oocyst counts from control group children

// input data
data {
  int<lower=0> N;
  array[N] int oocyst_count;
}

// parameters accepted by the model
parameters {
  real<lower=0> mu; // mean of neg bin distribution (untruncated, so includes all the 0s)
  real<lower=0> phi; // negative binomial dispersion par
}

// model to be estimated
model {
  #oocyst_count ~ neg_binomial_2(mu, phi);
  for (i in 1:N){
    real log_norm_const = log1m(neg_binomial_2_cdf( 0 | mu, phi )); // prob of x>0, on log scale; normalising constant

    if (oocyst_count[i] < 99) {
      target += neg_binomial_2_lpmf(oocyst_count[i] | mu, phi) - log_norm_const; // if oocysts count is <99 get the exact value density
    } else {
      target += neg_binomial_2_lccdf(98 | mu, phi) - log_norm_const; // prob of being at least 98
    }
  }

  mu ~ normal(10, 20); // approx because of mean and sd (oo_pos$oocyst_count)
  phi ~ exponential(2); // unsure, just basic exp
}

//
generated quantities{
    // simulations for PPC (removing the zeros and re-drawing)
    array[N] int sim_oocyst_counts;
    for ( i in 1:N ) sim_oocyst_counts[i] = neg_binomial_2_rng(mu, phi); // not truncated, will need to truncate later when comparing with data

    // log likelihood for comparison to other models
    vector[N] log_lik;
    // for ( i in 1:N ) log_lik[i] = neg_binomial_2_lpmf( oocyst_count[i] | mu, phi ) -
    //                               log1m(neg_binomial_2_cdf( 0 | mu, phi ));
    // subtracting log(1 - P(X=0)) from the neg bin likelihood to normalize it based on zero-trncation (same as 1/1-dnbinom(0,mu, phi))
    // and subtracting log(1 - P(x >98))
    for (i in 1:N) {
       real log_norm_const = log1m(neg_binomial_2_cdf(0 | mu, phi)); // same as above
       if (oocyst_count[i] < 99) {
         log_lik[i] = neg_binomial_2_lpmf(oocyst_count[i] | mu, phi) - log_norm_const;
       } else {
         log_lik[i] = neg_binomial_2_lccdf(98 | mu, phi) - log_norm_const;
       }
  }

}
// truncation for log-lik computes the log-density of the truncated distribution evaluated at the observed data points
