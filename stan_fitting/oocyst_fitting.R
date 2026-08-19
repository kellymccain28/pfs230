# Script to fit distribution of oocyst counts in DSF data from Patrick
library(bayesplot)
library(tidyverse)
library(rstan)
library(shinystan)

parallel::detectCores()
options(mc.cores = 4)

oo <- readxl::read_xlsx("C:/Users/kem22/OneDrive - Imperial College London/Challenger, Joseph's files - PfVIMT_Joe_Kelly/Mosquito_centric/AgeDE.Pfs230D1.DSF.Mosquitoes.xlsx")
oo <- oo %>%
  janitor::clean_names()

hist(oo[oo$oocyst_count > 0,]$oocyst_count)
ggplot(oo %>% filter(oocyst_count > 0)) +
  geom_histogram(aes(x = oocyst_count)) +
  theme_bw()

oo_pos <- oo %>%
  filter(oocyst_count > 0) %>%
  filter(vaccine_group == 'Comparator')
oocysts <- oo_pos$oocyst_count

oodatlist <- list(oocyst_count = oo_pos$oocyst_count,
                  N = nrow(oo_pos))

# First, run a model just estimating the parameters for a neg bin distribution that matches these data,
# removing the 0s in the model
# goal is to produce a PDF and a CDF to compare to teh one Joe had ade from the Bompard data

fit_control = rstan::stan('stan_fitting/stan_oocyst_dsf_control_children.stan',
                                model_name = 'control_nocov',
                         data = oodatlist,
                         iter = 5000,
                         chains = 4)
p <- c('mu', "phi")

fit_control_sml <- as.shinystan(fit_control, pars = p)
my_sso <- launch_shinystan(fit_control_sml)

print(fit_control, p)
mcmc_dens_overlay(fit_control, p)
mcmc_intervals(fit_control, p)

bayesplot::mcmc_trace(fit_control, p, nrow = 3)
# bayesplot::mcmc_trace(fit_control, pars = c("mu", "phi"), window = c(1000,1500), nrow = 3)

control_pairs <- bayesplot::mcmc_pairs(fit_control, p) # there is high correlation between the two parameters

np <- nuts_params(fit_control)
# mcmc_parcoord(fit_control, np = np, p)

samples_control <- rstan::extract(fit_control)# get samples from the posterior distribution
# dim(samples_control$mu)
ppccontrol <- bayesplot::ppc_dens_overlay(y = oo_pos$oocyst_count, yrep = samples_control$sim_oocyst_counts[1:100,]) +
  xlim(c(0,3))

# rhat
rhats <- rhat(fit_control, p)
mcmc_rhat(rhats, p)

# eff sample size
ratios <- neff_ratio(fit_control, p) # worry if ratio of neff to N is < 1
mcmc_neff(ratios)

# autocorrelation
bayesplot::mcmc_acf(fit_control, p) # slowly drops to 0 which indicates some autocorrelation


# Make CDF and PDF for fitted distribution
muhat <- mean(samples_control$mu)
phihat <- mean(samples_control$phi)

df_NB <- data.frame('o' = seq(1,100),
                    'pmf' = (1/(1-dnbinom(0, mu = muhat, size = phihat)))*dnbinom(1:100, mu = muhat, size = phihat),
                    'pmf_bompard' = (1/(1-dnbinom(0, mu = 0.000157, size = 0.00000495)))*dnbinom(1:100, mu = 0.000157, size = 0.00000495))
                    # above renormalizes the distribution after removing the 0s
df_NB$cmf <- cumsum(df_NB$pmf)
df_NB$cmf_bompard <- cumsum(df_NB$pmf_bompard)

pl0 <- ggplot(df_NB) +
  geom_point(aes(x = o, y = pmf, color = 'Fitted')) +
  geom_point(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
  scale_x_log10() + theme_classic() +
  #scale_y_log10() +
  xlab('Oocyst count') + ylab('Probability')+
  theme(legend.position = c(0.8, 0.7))
pl0

pl <- ggplot(df_NB) +
  geom_point(aes(x = o, y = pmf, color = 'Fitted')) +
  geom_point(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(x = 'Oocyst count',
       y = 'PMF')+
  theme(legend.position = c(0.3, 0.3))
pl

plc <- ggplot(df_NB) +
  geom_point(aes(x = o, y = cmf, colour = 'Fitted')) +
  geom_point(aes(x = o, y = cmf_bompard, colour = 'Bompard')) +
  scale_x_log10() +
  theme_classic() +
  scale_y_log10() +
  ylim(c(0,1)) +
  labs(x = 'Oocyst count',
       y ='Cumulative probability') +
  theme(legend.position = c(0.8, 0.3))
plc
cowplot::plot_grid(pl0, plc + theme(legend.position = 'none'))
ggsave('DSF_fitted_pmf.pdf', height = 5.2, width = 9)

############################################################################################################

# then add the vaccination group information into the data and update the model with a binary var for vax or not
# common effects / random effects? can try each of these and compare
ooall <- oo %>%
  filter(oocyst_count > 0)
ooall$vaccinated <- ifelse(ooall$vaccine_group == 'Comparator', 0, 1)
oo2datlist <- list(oocyst_count = ooall$oocyst_count,
                   vaccinated = ooall$vaccinated,
                  N = nrow(ooall))

fit_both = rstan::stan('stan_fitting/stan_oocyst_dsf_control_vaccinated_children.stan',
                          model_name = 'control_nocov',
                          data = oo2datlist,
                          iter = 5000,
                          chains = 4)
p <- c('mu0', 'beta_vacc', "phi")

fit_both_sml <- as.shinystan(fit_both, pars = p)
my_sso <- launch_shinystan(fit_both_sml)

print(fit_both, p)
mcmc_dens_overlay(fit_both, p)
mcmc_intervals(fit_both, p)

bayesplot::mcmc_trace(fit_both, p, nrow = 3)
# bayesplot::mcmc_trace(fit_both, pars = c("mu", "phi"), window = c(1000,1500), nrow = 3)

bothpairs <- bayesplot::mcmc_pairs(fit_both, p)

np <- nuts_params(fit_both)
# mcmc_parcoord(fit_both, np = np, p)

samples_both <- rstan::extract(fit_both)# get samples from the posterior distribution
# dim(samples_control$mu)
ppcboth <- bayesplot::ppc_dens_overlay(y = ooall$oocyst_count, yrep = samples_both$sim_oocyst_counts[1:100,]) +
  xlim(c(0,1))

# rhat
rhats <- rhat(fit_both, p)
mcmc_rhat(rhats, p)

# eff sample size
ratios <- neff_ratio(fit_both, p) # worry if ratio of neff to N is < 1
mcmc_neff(ratios)

# autocorrelation
bayesplot::mcmc_acf(fit_both, p) # slowly drops to 0 which indicates some autocorrelation

# Make CDF and PDF for fitted distribution -- neds modificaiton to do it for vaccinated and unvaccinated
muhat <- mean(samples_both$mu)
phihat <- mean(samples_both$phi)

df_NB <- data.frame('o' = seq(1,100),
                    'pmf' = (1/(1-dnbinom(0, mu = muhat, size = phihat)))*dnbinom(1:100, mu = muhat, size = phihat),
                    'pmf_bompard' = (1/(1-dnbinom(0, mu = 0.000157, size = 0.00000495)))*dnbinom(1:100, mu = 0.000157, size = 0.00000495))
# above renormalizes the distribution after removing the 0s
df_NB$cmf <- cumsum(df_NB$pmf)
df_NB$cmf_bompard <- cumsum(df_NB$pmf_bompard)

pl0 <- ggplot(df_NB) +
  geom_point(aes(x = o, y = pmf, color = 'Fitted')) +
  geom_point(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
  scale_x_log10() + theme_classic() +
  #scale_y_log10() +
  xlab('Oocyst count') + ylab('Probability')+
  theme(legend.position = c(0.8, 0.7))
pl0

pl <- ggplot(df_NB) +
  geom_point(aes(x = o, y = pmf, color = 'Fitted')) +
  geom_point(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(x = 'Oocyst count',
       y = 'PMF')+
  theme(legend.position = c(0.3, 0.3))
pl

plc <- ggplot(df_NB) +
  geom_point(aes(x = o, y = cmf, colour = 'Fitted')) +
  geom_point(aes(x = o, y = cmf_bompard, colour = 'Bompard')) +
  scale_x_log10() +
  theme_classic() +
  scale_y_log10() +
  ylim(c(0,1)) +
  labs(x = 'Oocyst count',
       y ='Cumulative probability') +
  theme(legend.position = c(0.8, 0.3))
plc
cowplot::plot_grid(pl0, plc + theme(legend.position = 'none'))
ggsave('DSF_fitted_pmf_vacc_control.pdf', height = 5.2, width = 9)
