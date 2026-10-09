# Script to run different scenarios assuming different models for the relationship between oocysts --> P(inf)
# which lead to different relationships between TRA-TBA+

library(malariasimulation)
library(postie)
library(ggplot2)

params <- get_parameters(overrides = list(
  human_population = 1000
))

params <- set_tbv(params,
                  timesteps = 365 * 12,
                  coverages = 0.5,
                  ages = 5:100,
                  adult_scaling = 1)

# Bompard parameters for oocyst distribtuion
mALL <- 0.000157   # NB mean
rrALL <- 0.00000495 # NB dispersion (size)

params$tbv_mu <- mALL
params$tbv_mt <- mALL
params$tbv_md <- mALL
params$tbv_ma <- mALL
params$tbv_k <- rrALL

sim <- run_simulation(timesteps = 365 * 15,
                      parameters = params)

sim <- drop_burnin(sim, burnin = 10*365)
proc <- get_prevalence(sim)

ggplot(proc) +
  geom_line(aes(x = time, y = lm_prevalence_2_10))+
  geom_smooth(aes(x = time, y = lm_prevalence_2_10))+
  geom_vline(aes(xintercept = 2000+2))
