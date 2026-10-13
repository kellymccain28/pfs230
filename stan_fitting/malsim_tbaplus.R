# Script to run different scenarios assuming different models for the relationship between oocysts --> P(inf)
# which lead to different relationships between TRA-TBA+

library(malariasimulation)
library(postie)
library(ggplot2)
library(tidyverse)

kmtheme <- theme_bw(base_size = 12)

# Bompard parameters for oocyst distribtuion
mALL <- 0.000157   # NB mean
rrALL <- 0.00000495 # NB dispersion (size)
year = 365

run_wtbaparams <- function(m, k,
                           tbv_timestep = 12*year,
                           tbv_coverage = 0.9,
                           tbv_ages =5:100,
                           timesteps = 20*year,
                           burnin = 10*year
                           ){

  if(timesteps < tbv_timestep){
    stop('TBV should be delivered before max timestep')}
  if(burnin > tbv_timestep){
    stop('TBV should be delivered after burnin')}

  params <- get_parameters(overrides = list(
    human_population = 10000
  )) %>%
    set_drugs(list(AL_params)) %>%
    set_clinical_treatment(1, 0.5, 0.5) %>%
    set_epi_outputs(
      age_group = c(0, 5, 15, 200) * year,
      prevalence = c(2, 10) * year,
      clinical_incidence = c(0, 5, 15, 200) * year,
      severe_incidence = c(0, 5, 15, 200) * year
    )

  params <- set_tbv(params,
                    timesteps = tbv_timestep,
                    coverages = tbv_coverage,
                    ages = tbv_ages,
                    adult_scaling = 1)

  params$tbv_mu <- m
  params$tbv_mt <- m
  params$tbv_md <- m
  params$tbv_ma <- m
  params$tbv_k <- k

  params <- set_equilibrium(params, init_EIR = 10)

  sim <- run_simulation(timesteps = timesteps,
                        parameters = params)
  sim <- drop_burnin(sim, burnin = burnin)

  proc <- get_prevalence(sim) %>%
    summarise(
      lm_prevalence_2_10 = mean(lm_prevalence_2_10),
      time = mean(time),
      .by = c(year, month)
    )
  inc <- get_rates(sim) %>%
    summarise(
      clinical = weighted.mean(clinical, person_days) * 365,
      severe = weighted.mean(severe, person_days) * 365,
      time = mean(time),
      .by = c(age_lower, age_upper, year, month)
    )

  infectivity_monthly <- sim %>%
    postie:::format_time(baseline_year = 2000) %>%
    select(starts_with('infectiv'), year, month) %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), mean),
      .by = c('month', 'year')
    )

  out <- left_join(inc, proc) %>%
    left_join(infectivity_monthly, by = c('year','month'))

  return(out)
}


proc_orig <- run_wtbaparams(m = mALL, k = rrALL) %>%
  mutate(model = 'm0')
proc_m1 <- run_wtbaparams(m = 0.3587039, k = 0.02142727) %>% # from tbaplus_fitting.R
  mutate(model = 'm1')
proc_m2 <- run_wtbaparams(m = 0.0001570028, k = 2.049708e-05) %>% # from tbaplus_fitting.R
  mutate(model = 'm2')

proc <- rbind(proc_orig, proc_m1, proc_m2)

cases_averted <- proc %>%
  select(age_lower, age_upper, clinical, model, time) %>%
  mutate(across(c(age_lower, age_upper),
                ~ round(.x, digits = 0))) %>%
  pivot_wider(names_from = c('model'),
              values_from = 'clinical') %>%
  mutate(ca_m1 = m0 - m1,
         ca_m2 = m0 - m2)

# Plot prevalence 2-10 per person per year, averaged by month
ggplot(proc) +
  geom_line(aes(x = time, y = lm_prevalence_2_10, color = model)) +
  geom_vline(aes(xintercept = 2000+2), color = 'red', linetype = 2) +
  # facet_wrap(~model) +
  kmtheme


# Plot monthly incidence by age group
ggplot(proc) +
  geom_line(aes(x = time, y = clinical, color = model)) +
  facet_wrap(~age_lower, nrow = 3) +
  kmtheme +
  theme(legend.position = 'bottom')

# Plot incidence averted by age group
ggplot(cases_averted %>% filter(year < 2005)) +
  geom_line(aes(x = time, y = ca_m1, color = 'm1')) +
  geom_line(aes(x = time, y = ca_m2, color = 'm2')) +
  geom_vline(aes(xintercept = 2000+2), color = 'red', linetype = 2) +
  geom_hline(aes(yintercept = 0), color = 'darkgreen', linetype = 2) +
  labs(y = 'cases averted') +
  facet_wrap(~age_lower,
             nrow = 3) +
  kmtheme

# Plot infectivity
ggplot(proc) +
  geom_line(aes(x = time, y = infectivity, color = model)) +
  kmtheme +
  theme(legend.position = 'bottom')
