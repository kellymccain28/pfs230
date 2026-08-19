# Script to fit a model to step 1 (oocysts -> total sporozoites)
library(rstan)
parallel::detectCores()
options(mc.cores = 4)

d1_spz_oocysts <- readxl::read_xlsx("C:/Users/kem22/OneDrive - Imperial College London/Challenger, Joseph's files - PfVIMT_Joe_Kelly/literature-oocysts-spz.xlsx",
                                               sheet = 3) %>%
  mutate(oocysts = as.integer(oocysts),
         spz_total = as.integer(spz_total)) %>%
  filter(oocysts > 0)


#Fit saturating stan model ----
dat_list <- list(oocysts = d1_spz_oocysts$oocysts,
                 spz_total = d1_spz_oocysts$spz_total,
                 N = nrow(d1_spz_oocysts))
fit_sat = rstan::stan('stan_fitting/stan_step1_saturating.stan',
                  data = dat_list,
                  iter = 3000,
                  chains = 4)
saveRDS(fit_sat, 'stan_fitting/step1_fit_sat.rds')
print(fit_sat, pars = c("alpha", "beta", "phi", "mu[1]"))
plot(fit_sat, pars = c("alpha", "beta", "phi"))
rstan::traceplot(fit_sat, pars = c("alpha", "beta", "phi"), nrow = 3)
pairssat <- bayesplot::mcmc_pairs(fit_sat, pars = c("alpha", "beta", "phi", "mu[1]"))
ggsave('stan_fitting/step1_fit_sat_pairs.png', pairssat, height = 6, width = 6)

bayesplot::mcmc_acf(fit_sat, p = c("alpha", "beta", "phi"))

samples_sat <- rstan::extract(fit_sat)# get samples from the posterior distribution
dim(samples_sat$mu)

spz_bar <- mean(d1_spz_oocysts$spz_total) #

# fitted mean
ggplot(d1_spz_oocysts) +
  geom_line(aes(x=oocysts, y=colMeans(samples_sat$mu),
             color= 'predicted'), linewidth = 1) +
  geom_point(aes(x=oocysts, y=spz_total,
             color= 'observed'), size = 2) +
  labs(x = 'oocysts',
       y = 'total spz') +
  scale_y_log10() +
  # scale_x_log10() +
  theme_classic()
ggsave('stan_fitting/step1_fit_sat_fittedmean.png', height = 6, width = 6)

# to do
# ppc checks, with bayesplot
ppcsat <- bayesplot::ppc_dens_overlay(y = d1_spz_oocysts$spz_total, yrep = samples_sat$sim_values_tot_spz[1:100,])
ggsave('stan_fitting/step1_fit_sat_ppc.png', ppcsat, height = 6, width = 6)


###############################################################################################

#Fit flat stan model ----
dat_list <- list(oocysts = d1_spz_oocysts$oocysts,
                 spz_total = d1_spz_oocysts$spz_total,
                 N = nrow(d1_spz_oocysts))
fit_flat = rstan::stan('stan_fitting/stan_step1_flat.stan',
                  data = dat_list,
                  iter = 3000,
                  chains = 4)
saveRDS(fit_flat, 'stan_fitting/step1_fit_flat.rds')
print(fit_flat, pars = c('mu', "phi"))
plot(fit_flat, pars = c("mu", "phi"))
rstan::traceplot(fit_flat, pars = c("mu", "phi"), nrow = 3)
pairsflat <- bayesplot::mcmc_pairs(fit_flat, pars = c("mu", "phi"))
ggsave('stan_fitting/step1_fit_flat_pairs.png', pairsflat, height = 6, width = 6)

bayesplot::mcmc_acf(fit_flat, p = c("mu", "phi"))

samples_flat <- rstan::extract(fit_flat)# get samples from the posterior distribution
dim(samples_flat$mu)

spz_bar <- mean(d1_spz_oocysts$spz_total) #

# fitted mean
ggplot(d1_spz_oocysts) +
  geom_line(aes(x=oocysts, y=mean(samples_flat$mu),
                color= 'predicted'), linewidth = 1) +
  geom_point(aes(x=oocysts, y=spz_total,
                 color= 'observed'), size = 2) +
  labs(x = 'oocysts',
       y = 'total spz') +
  scale_y_log10() +
  # scale_x_log10() +
  theme_classic()
ggsave('stan_fitting/step1_fit_flat_fittedmean.png', height = 6, width = 6)

# to do
# ppc checks, with bayesplot
ppcflat <- bayesplot::ppc_dens_overlay(y = d1_spz_oocysts$spz_total, yrep = samples_flat$sim_values_tot_spz[1:100,])
ggsave('stan_fitting/step1_fit_flat_ppc.png', ppcflat, height = 6, width = 6)

###############################################################################################

#Fit linear stan model ----
dat_list <- list(oocysts = d1_spz_oocysts$oocysts,
                 spz_total = d1_spz_oocysts$spz_total,
                 N = nrow(d1_spz_oocysts))
fit_linear = rstan::stan('stan_fitting/stan_step1_linear.stan',
                       data = dat_list,
                       iter = 3000,
                       chains = 4)
saveRDS(fit_linear, 'stan_fitting/step1_fit_linear.rds')

print(fit_linear, pars = c('m', "phi", "mu[1]"))
plot(fit_linear, pars = c("m", "phi"))
rstan::traceplot(fit_linear, pars = c("m", "phi"), nrow = 3)
pairslin <- bayesplot::mcmc_pairs(fit_linear, pars = c("m", "phi", "mu[1]"))
ggsave('stan_fitting/step1_fit_linear_pairs.png', pairslin, height = 6, width = 6)

bayesplot::mcmc_acf(fit_linear, p = c("m", "phi"))

samples_linear <- rstan::extract(fit_linear)# get samples from the posterior distribution
dim(samples_linear$mu)

spz_bar <- mean(d1_spz_oocysts$spz_total) #

# fitted mean
ggplot(d1_spz_oocysts) +
  geom_line(aes(x=oocysts, y=colMeans(samples_linear$mu),
                color= 'predicted'), linewidth = 1) +
  geom_point(aes(x=oocysts, y=spz_total,
                 color= 'observed'), size = 2) +
  labs(x = 'oocysts',
       y = 'total spz') +
  scale_y_log10() +
  # scale_x_log10() +
  theme_classic()
ggsave('stan_fitting/step1_fit_linear_fittedmean.png', height = 6, width = 6)

# to do
# ppc checks, with bayesplot
ppclinear <- bayesplot::ppc_dens_overlay(y = d1_spz_oocysts$spz_total, yrep = samples_linear$sim_values_tot_spz[1:100,])
ggsave('stan_fitting/step1_fit_linear_ppc.png', ppclinear, height = 6, width = 6)


# Compare models
library(loo)

ll_sat <- loo::loo(loo::extract_log_lik(fit_sat))
ll_flat <- loo::loo(loo::extract_log_lik(fit_flat))
ll_linear <- loo::loo(loo::extract_log_lik(fit_linear))

loo::loo_compare(list(saturated = ll_sat, flat = ll_flat, linear = ll_linear))

rethinking::compare(fit_sat,fit_flat,fit_linear, func = 'WAIC')
