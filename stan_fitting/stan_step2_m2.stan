// Stan model to fit the step 2 of Joe's model
// (ii) Total SPZ -> SPZ per bite
// g3 <- function(x) 50 + 4950 * (x**2/(x**2 + 10000**2))

// input data
data {
  int<lower=0> N;
  array[N] int spz_total;
  array[N] real spz_per_bite;
}

// parameters accepted by the model
parameters {
  real<lower=0> sigma; // sd of underlying normal distribution of log(spz_per_bite)
  real<lower=0> a; //
  real<lower=0> b; //
  real<lower=0> c; //
}

// transformed parameters
// transformed parameters {
//   vector[N] mu = a + b * (to_vector(spz_total) ^2 ./ ( to_vector(spz_total)^2 + c^2)); // mean of underlying normal distribution of log(spz_per_bite)
//   // above uses vectorized form which is integers; vector allows for computations
// }

// model to be estimated
model {
  vector[N] mu = a + b * (to_vector(spz_total) ^2 ./ ( to_vector(spz_total)^2 + c^2)); // mean of underlying normal distribution of log(spz_per_bite)
  spz_per_bite ~ lognormal(log(mu), sigma);

  a ~ normal(50, 10); //
  b ~ normal(5000, 100);
  c ~ normal(10000, 1000);
  sigma ~ exponential(0.01);
}

//
generated quantities{
    vector[N] mu = a + b * (to_vector(spz_total) ^2 ./ ( to_vector(spz_total)^2 + c^2)); // mean of underlying normal distribution of log(spz_per_bite)

    // simulations for PPC
    array[N] real sim_spz_per_bite;
    for ( i in 1:N ) sim_spz_per_bite[i] = lognormal_rng(log(mu[i]), sigma);

    // log likelihood for comparison to other models
    vector[N] log_lik;
    for ( i in 1:N ) log_lik[i] = lognormal_lpdf( spz_per_bite[i] | log(mu[i]), sigma );


}
