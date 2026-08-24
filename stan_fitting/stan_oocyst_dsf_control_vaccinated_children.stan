// Model of oocyst counts comparing control and vaccinated children
// where vaccination is represented with a binary variable

// input data
data {
  int<lower=0> N;
  array[N] int oocyst_count;
  array[N] int<lower=0, upper=1> vaccinated; // binary variable for vaccinated or unvaccinated
}

// parameters accepted by the model
parameters {
  real<lower=0> mu0; // baseline mean of neg bin distribution
  real beta_vacc; // vaccine effect
  real<lower=0> phi; // negative binomial dispersion par
}

// transformed parameters
transformed parameters {
  vector[N] mu;
  for ( i in 1:N ) mu[i] = mu0 * exp(beta_vacc * vaccinated[i]); // mean of neg bin distribution
  // above uses vectorized form of the oocysts assay which is integers; vector allows for computations
}

// model to be estimated
model {
  // oocyst_count ~ neg_binomial_2(mu, phi);# T[1 , 99];
  for (i in 1:N) {
    real log_norm_const = log1m(neg_binomial_2_cdf(0 | mu[i], phi)); // P(X > 0), on log scale

    if (oocyst_count[i] < 99) {
      target += neg_binomial_2_lpmf(oocyst_count[i] | mu[i], phi) - log_norm_const;
    } else {
      target += neg_binomial_2_lccdf(98 | mu[i], phi) - log_norm_const;
    }
  }

  mu0 ~ normal(10, 20); // approx because of mean and sd (oo_pos$oocyst_count)
  phi ~ exponential(0.01); // unsure, just basic exp
  beta_vacc ~ normal(0, 1); // centered around 0 which means that there is no difference; neg would be reduction in oocysts w/ vaccination
}

//
generated quantities{
    // simulations for PPC (removing the zeros and re-drawing)
    array[N] int sim_oocyst_counts;
    for ( i in 1:N ) sim_oocyst_counts[i] = neg_binomial_2_rng(mu[i], phi);

    // log likelihood for comparison to other models
    // truncation for log-lik computes the log-density of the truncated distribution evaluated at the observed data points
    vector[N] log_lik;
    // for ( i in 1:N ) log_lik[i] = neg_binomial_2_lpmf( oocyst_count[i] | mu[i], phi ) -
    //                               log1m(neg_binomial_2_cdf( 0 | mu[i], phi ));
    for (i in 1:N) {
      real log_norm_const = log1m(neg_binomial_2_cdf(0 | mu[i], phi));

      if (oocyst_count[i] < 99) {
        log_lik[i] = neg_binomial_2_lpmf(oocyst_count[i] | mu[i], phi) - log_norm_const;
      } else {
        log_lik[i] = neg_binomial_2_lccdf(98 | mu[i], phi) - log_norm_const;
      }
    }
  // subtracting log(1 - P(X=0)) from the neg bin likelihood to normalize it based on zero-trncation (same as 1/1-dnbinom(0,mu, phi))
    // will also need to treat the oocyst counts with 99 or above differently becuase they were reocrded as at least 99 (not truncated but censored)?

}

