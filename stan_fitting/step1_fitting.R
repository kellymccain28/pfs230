# Script to fit a model to step 1 (oocysts -> total sporozoites)
library(rstan)

d1_spz_oocysts <- readxl::read_xlsx("C:/Users/kem22/OneDrive - Imperial College London/Challenger, Joseph's files - PfVIMT_Joe_Kelly/literature-oocysts-spz.xlsx",
                                               sheet = 3) %>%
  mutate(oocysts = as.integer(oocysts),
         spz_total = as.integer(spz_total)) %>%
  filter(oocysts > 0)


#Fit stan model
parallel::detectCores()
options(mc.cores = 4)

dat_list <- list(oocysts = d1_spz_oocysts$oocysts,
                 spz_total = d1_spz_oocysts$spz_total,
                 N = nrow(d1_spz_oocysts))
fit = rstan::stan('stan_fitting/stan_step1.stan',
                  data = dat_list,
                  iter = 3000,
                  chains = 4)
print(fit, pars = c("alpha", "beta", "phi"))
plot(fit, pars = c("alpha", "beta", "phi"))
rstan::traceplot(fit, pars = c("alpha", "beta", "phi"), nrow = 3)
pairs(fit, pars = c("alpha", "beta", "phi"))

samples <- rstan::extract(fit)# get samples from the posterior distribution
dim(samples$mu)

spz_bar <- mean(d1_spz_oocysts$spz_total) #

# below is using the average value of beta and alpha but this is not corret
# mean(samples$alpha)
# mean(samples$beta)
# # f <- function(oocysts) (oocysts * mean(samples$beta)) / (1 + oocysts / mean(samples$alpha))
# ggplot(d1_spz_oocysts) +
#   geom_line(aes(x = oocysts, y = f(oocysts))) +
#   scale_y_log10() +
#   theme_classic()

# fitted mean
ggplot(d1_spz_oocysts) +
  geom_line(aes(x=oocysts, y=colMeans(samples$mu),
             color= 'predicted')) +
  geom_point(aes(x=oocysts, y=spz_total,
             color= 'observed')) +
  scale_y_log10() +
  theme_classic()

# to do
# ppc checks, with bayesplot
bayesplot::ppc_dens_overlay(y = d1_spz_oocysts$spz_total, yrep = samples$sim_values_tot_spz)
