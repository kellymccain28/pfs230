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
saveRDS(fit_control, 'stan_fitting/outputs/oocyst_control.rds')
p <- c('mu', "phi")

# fit_control_sml <- as.shinystan(fit_control, pars = p)
# my_sso <- launch_shinystan(fit_control_sml)

print(fit_control, p)
mcmc_dens_overlay(fit_control, p)
mcmc_intervals(fit_control, p)

bayesplot::mcmc_trace(fit_control, p, nrow = 3)
# bayesplot::mcmc_trace(fit_control, pars = c("mu", "phi"), window = c(1000,1500), nrow = 3)

control_pairs <- bayesplot::mcmc_pairs(fit_control, p) # there is high correlation between the two parameters
ggsave('stan_fitting/outputs/oocysts_fit_control_pairs.png', control_pairs, height = 6, width = 6)

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


# Plot fitted mean compared to the data
# first, need to get the mean of the fitted distribution, but truncated
sim_oocysts_matrix <- samples_control$sim_oocyst_counts # dimension: draws x N where draws = iter/2 * nchains
mean_fitted <- rowSums(sim_oocysts_matrix * (sim_oocysts_matrix > 0)) / rowSums(sim_oocysts_matrix > 0)

# how many draws are NaN per group, worth reporting alongside the plot
n_na <- sum(is.nan(mean_fitted))
n_na
# summarise each group's truncated fitted mean (posterior mean + credible interval)
trunc_summary <- data.frame(
  vaccine_group = c('Comparator'),
  mu_fit_trunc = mean(mean_fitted, na.rm = TRUE),
  lower = quantile(mean_fitted, 0.025, na.rm = TRUE),
  upper = quantile(mean_fitted, 0.975, na.rm = TRUE)
)

controlfitted <- ggplot(oo_pos) +
  geom_jitter(aes(x = vaccine_group, y = oocyst_count, color = vaccine_group),
              width = 0.2, alpha = 0.5) +
  geom_hline(aes(yintercept = mean(oocyst_count)), linetype = 2, color = 'grey50') +
  geom_pointrange(data = trunc_summary,
                  aes(x = vaccine_group, y = mu_fit_trunc, ymin = lower, ymax = upper),
                  color = 'black', size = 0.8) +
  scale_y_log10() +
  theme_classic(base_size = 12) +
  labs(y = 'Oocyst count', x = 'Vaccine group',
       color = 'Vaccine group',
       caption = paste0('Truncated mean excludes draws with all-zero simulated counts\n(',
                        'Comparator: ', n_na0, '/', nrow(truncated_means), ' draws excluded)')) +
  theme(legend.position = c(0.8,0.8))
ggsave('stan_fitting/outputs/oocysts_fit_control_fittedmeans.png', controlfitted, height = 6, width = 6)


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
  geom_line(aes(x = o, y = pmf, color = 'Fitted')) +
  geom_line(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
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
  geom_line(aes(x = o, y = pmf, color = 'Fitted')) +
  geom_line(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
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
  geom_line(aes(x = o, y = cmf, colour = 'Fitted')) +
  geom_line(aes(x = o, y = cmf_bompard, colour = 'Bompard')) +
  scale_x_log10() +
  theme_classic() +
  scale_y_log10() +
  ylim(c(0,1)) +
  labs(x = 'Oocyst count',
       y ='Cumulative probability') +
  theme(legend.position = c(0.8, 0.3))
plc
cowplot::plot_grid(pl0, plc + theme(legend.position = 'none'))
ggsave('stan_fitting/outputs/DSF_fitted_pmf.png', height = 5.2, width = 9)

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
saveRDS(fit_both, 'stan_fitting/outputs/oocyst_control_vaccinated.rds')
p <- c('mu0', 'beta_vacc', "phi")


# fit_both_sml <- as.shinystan(fit_both, pars = p)
# my_sso <- launch_shinystan(fit_both_sml)

print(fit_both, p)
# sense check that the prop of zeroes is reasonably similar to data
oo %>% group_by(vaccine_group) %>% summarize(mean(oocyst_count == 0))
dnbinom(0, mu = 0.38, size = 0.01) # unvaccinated
dnbinom(0, mu = 0.38 * exp(-2.92), size = 0.01) # vaccinated

mcmc_dens_overlay(fit_both, p)
mcmc_intervals(fit_both, p)

bayesplot::mcmc_trace(fit_both, p, nrow = 3)
# bayesplot::mcmc_trace(fit_both, pars = p, window = c(1000,1200), nrow = 3)

bothpairs <- bayesplot::mcmc_pairs(fit_both, p)
ggsave('stan_fitting/outputs/oocysts_fit_both_pairs.png', bothpairs, height = 6, width = 6)

np <- nuts_params(fit_both)
mcmc_parcoord(fit_both, np = np, p)

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

# divergence
mcmc_nuts_divergence(np, )

# autocorrelation
acfboth <- bayesplot::mcmc_acf(fit_both, p) # slowly drops to 0 which indicates some autocorrelation
ggsave('stan_fitting/outputs/oocysts_fit_both_autocorrelation.png', acfboth, height = 6, width = 6)


# Plot fitted mean compared to the data
# first, need to get the mean of the fitted distribution, but truncated
sim_oocysts_matrix <- samples_both$sim_oocyst_counts # dimension: draws x N where draws = iter/2 * nchains
# truncated mean per group by draw:
simcontrol <- sim_oocysts_matrix[, ooall$vaccinated == 0]
simvaccinated <- sim_oocysts_matrix[, ooall$vaccinated == 1]

mean_control <- rowSums(simcontrol * (simcontrol > 0)) / rowSums(simcontrol > 0)
mean_vaccinated <- rowSums(simvaccinated * (simvaccinated > 0)) / rowSums(simvaccinated > 0)

# how many draws are NaN per group, worth reporting alongside the plot
n_na0 <- sum(is.nan(mean_control))
n_na1 <- sum(is.nan(mean_vaccinated))
n_na0; n_na1
# summarise each group's truncated fitted mean (posterior mean + credible interval)
trunc_summary <- data.frame(
  vaccine_group = c('Comparator', 'Pfs230'),
  mu_fit_trunc = c(mean(mean_control, na.rm = TRUE), mean(mean_vaccinated, na.rm = TRUE)),
  lower = c(quantile(mean_control, 0.025, na.rm = TRUE), quantile(mean_vaccinated, 0.025, na.rm = TRUE)),
  upper = c(quantile(mean_control, 0.975, na.rm = TRUE), quantile(mean_vaccinated, 0.975, na.rm = TRUE)),
  mean_obs = c(mean(ooall[ooall$vaccinated==0,]$oocyst_count),mean(ooall[ooall$vaccinated==1,]$oocyst_count))
)

bothfitted <- ggplot(ooall) +
  geom_jitter(aes(x = vaccine_group, y = oocyst_count, color = vaccine_group),
              width = 0.2, alpha = 0.5) +
  geom_hline(data = trunc_summary,
             aes(yintercept = mean_obs, group = vaccine_group, color = vaccine_group), linetype = 2) +
  geom_pointrange(data = trunc_summary,
                  aes(x = vaccine_group, y = mu_fit_trunc, ymin = lower, ymax = upper),
                  color = 'black', size = 0.5) +
  scale_y_log10() +
  theme_classic(base_size = 12) +
  labs(y = 'Oocyst count', x = 'Vaccine group',
       color = 'Vaccine group',
       caption = paste0('Truncated mean excludes draws with all-zero simulated counts\n(',
                        'Comparator: ', n_na0, '/', nrow(truncated_means), ' draws excluded; ',
                        'Pfs230: ', n_na1, '/', nrow(truncated_means), ' draws excluded)')) +
  theme(legend.position = c(0.8,0.8))
ggsave('stan_fitting/outputs/oocysts_fit_both_fittedmeans.png', bothfitted, height = 6, width = 6)


# Make CDF and PDF for fitted distribution -- neds modificaiton to do it for vaccinated and unvaccinated
# need to calculate mu for the vaccinated group, then take mean of that; and can use mu0 mean for unvacc (but I already have mu saved directly)
muhat <- data.frame(vaccinated = ooall$vaccinated,
                    mu_fit = colMeans(samples_both$mu)) %>%
  distinct()
phihat <- mean(samples_both$phi)

df_NB <- data.frame('o' = seq(1,100),
                    'pmf_vacc' = (1/(1-dnbinom(0, mu = muhat[muhat$vaccinated==1,]$mu_fit, size = phihat)))*
                      dnbinom(1:100, mu = muhat[muhat$vaccinated==1,]$mu_fit, size = phihat),
                    'pmf_control' = (1/(1-dnbinom(0, mu = muhat[muhat$vaccinated==0,]$mu_fit, size = phihat)))*
                      dnbinom(1:100, mu = muhat[muhat$vaccinated==0,]$mu_fit, size = phihat),
                    'pmf_bompard' = (1/(1-dnbinom(0, mu = 0.000157, size = 0.00000495)))*dnbinom(1:100, mu = 0.000157, size = 0.00000495))
# above renormalizes the distribution after removing the 0s
df_NB$cmf_vacc <- cumsum(df_NB$pmf_vacc)
df_NB$cmf_control <- cumsum(df_NB$pmf_control)
df_NB$cmf_bompard <- cumsum(df_NB$pmf_bompard)

df_NB_long <- df_NB %>%
  pivot_longer(cols = !o,
               names_to = c('pr_type','group'),
               names_pattern = c('(.*)_(.*)'),
               values_to = 'value')
ggplot(df_NB_long %>% filter(pr_type == 'pmf')) +
  geom_col(aes(x = o, y = value, group = group, fill = group),
           position = position_dodge()) +
  facet_wrap(~group) + theme_classic()
ggplot(df_NB_long %>% filter(pr_type == 'cmf')) +
  geom_col(aes(x = o, y = value, group = group, fill = group),
           position = position_dodge()) +
  facet_wrap(~group) + theme_classic()

pl0 <- ggplot(df_NB) +
  geom_line(aes(x = o, y = pmf_vacc, color = 'Fitted vaccine group')) +
  geom_line(aes(x = o, y = pmf_control, color = 'Fitted control group')) +
  geom_line(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
  geom_point(aes(x = o, y = pmf_vacc, color = 'Fitted vaccine group')) +
  geom_point(aes(x = o, y = pmf_control, color = 'Fitted control group')) +
  geom_point(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
  scale_x_log10() + theme_classic() +
  #scale_y_log10() +
  xlab('Oocyst count') + ylab('Probability')+
  theme(legend.position = c(0.8, 0.7))
pl0

pl <- ggplot(df_NB) +
  geom_line(aes(x = o, y = pmf_vacc, color = 'Fitted vaccine group')) +
  geom_line(aes(x = o, y = pmf_control, color = 'Fitted control group')) +
  geom_line(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
  geom_point(aes(x = o, y = pmf_vacc, color = 'Fitted vaccine group')) +
  geom_point(aes(x = o, y = pmf_control, color = 'Fitted control group')) +
  geom_point(aes(x = o, y = pmf_bompard, color = 'Bompard')) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(x = 'Oocyst count',
       y = 'PMF')+
  theme(legend.position = c(0.3, 0.3))
pl

plc <- ggplot(df_NB) +
  geom_line(aes(x = o, y = cmf_vacc, colour = 'Fitted vaccine group')) +
  geom_line(aes(x = o, y = cmf_control, colour = 'Fitted control group')) +
  geom_line(aes(x = o, y = cmf_bompard, colour = 'Bompard')) +
  geom_point(aes(x = o, y = cmf_vacc, colour = 'Fitted vaccine group')) +
  geom_point(aes(x = o, y = cmf_control, colour = 'Fitted control group')) +
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
ggsave('stan_fitting/outputs/DSF_fitted_pmf_vacc_control.png', height = 5.2, width = 9)
