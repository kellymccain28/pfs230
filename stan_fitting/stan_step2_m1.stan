// Stan model to fit the step 2 of Joe's model
// (ii) Total SPZ -> SPZ per bite
// g6

// input data
data {
  int<lower=0> N;
  array[N] int spz_total;
  array[N] real spz_per_bite;
}

// parameters accepted by the model
parameters {
  real<lower=0> sigma; // sd of underlying normal distribution of log(spz_per_bite)
  real<lower=0, upper=1> floor_p; // floor prop of total spz ejected per bite
  real<lower=0> scale; // scale parameter
  real<lower=0, upper=1> mult; // multiplier
}

// transformed parameters
transformed parameters {
  vector[N] mu = to_vector(spz_total) .* ( floor_p + mult * 1 / ( 1 + to_vector(spz_total) / scale ) ); // mean of underlying normal distribution of log(spz_per_bite)
  // above uses vectorized form which is integers; vector allows for computations
}

// model to be estimated
model {
  spz_per_bite ~ lognormal(log(mu), sigma);

  floor_p ~ beta(1, 1.5); //
  scale ~ normal(40000, 10000);
  mult ~ beta(1, 1);
  sigma ~ exponential(0.01);
}

//
