# Script for Step 2 fitting
d2_dat <- readxl::read_xlsx("C:/Users/kem22/OneDrive - Imperial College London/Challenger, Joseph's files - PfVIMT_Joe_Kelly/literature-oocysts-spz.xlsx",
                                                     sheet = 4) %>%
  filter(mosquito == 'An coluzzi') %>%
  select(spz_total, spz_per_bite)


#Fit first stan model ---- g6
dat_list <- list(spz_total = d2_dat$spz_total,
                 spz_per_bite = d2_dat$spz_per_bite,
                 N = nrow(d2_dat))
fit_m1 = rstan::stan('stan_fitting/stan_step2_m1.stan',
                     data = dat_list,
                     iter = 3000,
                     chains = 4)
saveRDS(fit_m1, 'stan_fitting/step2_fit_m1.rds')

p = c('mu[1]','sigma','floor_p','scale','mult')

print(fit_m1, pars = p)
plot(fit_m1, pars = p)
rstan::traceplot(fit_m1, pars = p, nrow = 3)
pairsm1 <- bayesplot::mcmc_pairs(fit_m1, pars = p)
ggsave('stan_fitting/step2_fit_m1_pairs.png', pairsm1, height = 6, width = 6)

acfm1 <- bayesplot::mcmc_acf(fit_m1, p = p)
ggsave('stan_fitting/step2_fit_m1_acf.png', acfm1, height = 6, width = 6)

samples_m1 <- rstan::extract(fit_m1)# get samples from the posterior distribution

# fitted mean
ggplot(d2_dat) +
  geom_line(aes(x=spz_total, y=colMeans(samples_m1$mu),
                color= 'predicted'), linewidth = 1) +
  geom_point(aes(x=spz_total, y=spz_per_bite,
                 color= 'observed'), size = 2) +
  labs(x = 'total spz',
       y = 'spz per bite') +
  scale_y_log10() +
  # scale_x_log10() +
  theme_classic()
ggsave('stan_fitting/step2_fit_m1_fittedmean.png', height = 6, width = 6)

# to do
# ppc checks, with bayesplot
ppcm1 <- bayesplot::ppc_dens_overlay(y = log(d2_dat$spz_per_bite), yrep = log(samples_m1$sim_spz_per_bite[1:100,]))
ggsave('stan_fitting/step2_fit_m1_ppc.png', ppcm1, height = 6, width = 6)


#Fit second stan model ---- g3
dat_list <- list(spz_total = d2_dat$spz_total,
                 spz_per_bite = d2_dat$spz_per_bite,
                 N = nrow(d2_dat))
fit_m2 = rstan::stan('stan_fitting/stan_step2_m2.stan',
                     data = dat_list,
                     iter = 3000,
                     chains = 4)
saveRDS(fit_m2, 'stan_fitting/step2_fit_m2.rds')

p = c('mu[1]','sigma','a','b','c')

print(fit_m2, pars = p)
plot(fit_m2, pars = p)
rstan::traceplot(fit_m2, pars = p, nrow = 3)
pairsm2 <- bayesplot::mcmc_pairs(fit_m2, pars = p)
ggsave('stan_fitting/step2_fit_m2_pairs.png', pairsm2, height = 6, width = 6)

acfm2 <- bayesplot::mcmc_acf(fit_m2, p = p)
ggsave('stan_fitting/step2_fit_m2_acf.png', acfm2, height = 6, width = 6)

samples_m2 <- rstan::extract(fit_m2)# get samples from the posterior distribution

# fitted mean
ggplot(d2_dat) +
  geom_line(aes(x=spz_total, y=colMeans(samples_m2$mu),
                color= 'predicted'), linewidth = 1) +
  geom_point(aes(x=spz_total, y=spz_per_bite,
                 color= 'observed'), size = 2) +
  labs(x = 'total spz',
       y = 'spz per bite') +
  scale_y_log10() +
  # scale_x_log10() +
  theme_classic()
ggsave('stan_fitting/step2_fit_m2_fittedmean.png', height = 6, width = 6)

# to do
# ppc checks, with bayesplot
ppcm2 <- bayesplot::ppc_dens_overlay(y = log(d2_dat$spz_per_bite), yrep = log(samples_m2$sim_spz_per_bite[1:100,]))
ggsave('stan_fitting/step2_fit_m2_ppc.png', ppcm2, height = 6, width = 6)


##### Compare two models ----
loom1 <- loo::loo(loo::extract_log_lik(fit_m1))
loom2 <- loo::loo(loo::extract_log_lik(fit_m2))

loo_compare(list(m1 = loom1, m2 = loom2))

rethinking::compare(fit_m1,fit_m2, func = 'WAIC')
