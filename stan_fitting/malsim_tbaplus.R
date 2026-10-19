# Script to run different scenarios assuming different models for the relationship between oocysts --> P(inf)
# which lead to different relationships between TRA-TBA+
# pak::pak("mrc-ide/malariasimulation")
library(malariasimulation) # make sure this is the version on the main branch
library(postie)
library(ggplot2)
library(tidyverse)

model_cols <- c('#ff595e','#8ac926','#6a4c93')
kmtheme <- theme_bw(base_size = 12)
# Bompard parameters for oocyst distribtuion
mALL <- 0.000157   # NB mean
rrALL <- 0.00000495 # NB dispersion (size)
year = 365

m0_pars <- c(model = 'm0', m = mALL, k = rrALL)
m1_pars <- c(model = 'm1', m = 0.3587039, k = 0.02142727)
m2_pars <- c(model = 'm2', m = 0.0001570028, k = 2.049708e-05)
model_pars <- data.frame(rbind(m0_pars, m1_pars, m2_pars)) %>%
  mutate(m = as.numeric(m),
         k = as.numeric(k))
rownames(model_pars) <- NULL

run_wtbaparams <- function(m, k,
                           pop,
                           draw,
                           season = 'seasonal',
                           tbv_timestep = 12*year,
                           tbv_coverage = 0.9,
                           tbv_ages = 5:200,
                           timesteps = 20*year,
                           burnin = 10*year
                           ){

  if(timesteps < tbv_timestep){
    stop('TBV should be delivered before max timestep')}
  if(burnin > tbv_timestep){
    stop('TBV should be delivered after burnin')}

  seas_name <- 'seasonal'
  seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
  s2 <- tibble(seasonality, seas_name)

  seas_name <- 'perennial'
  seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
  s3 <- tibble(seasonality, seas_name)
  seasonality <- bind_rows(s2, s3) %>%
    filter(seas_name == season)

  params <- get_parameters(overrides = list(
    human_population = pop,
    model_seasonality = TRUE,
    # rainfall fourier parameters
    g0 = as.numeric(unlist(seasonality)[1]),
    g = as.numeric(unlist(seasonality)[2:4]),
    h = as.numeric(unlist(seasonality)[5:7]),
    individual_mosquitoes = FALSE
  )) %>%
    set_drugs(list(AL_params)) %>%
    set_clinical_treatment(1, 0.5, 0.5) %>%
    set_epi_outputs(
      age_group = c(0, 5, 15, 200) * year,
      prevalence = c(2, 10) * year,
      clinical_incidence = c(0, 5, 15, 200) * year,
      severe_incidence = c(0, 5, 15, 200) * year
    ) %>%
    set_parameter_draw(draw)

  peak <- peak_season_offset(params)
  tbv_timestep_offset <- tbv_timestep + peak - (4 * round(365/12))

  params <- set_tbv(params,
                    timesteps = tbv_timestep_offset, # 4 months before peak
                    coverages = tbv_coverage,
                    ages = tbv_ages)

  params$tbv_mu <- m
  params$tbv_mt <- m
  params$tbv_md <- m
  params$tbv_ma <- m
  params$tbv_k <- k

  params <- set_equilibrium(params, init_EIR = 10)

  set.seed(draw) # to allow for paired comparisons between scenarios
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
      clinical = weighted.mean(clinical, person_days) * 365/12,
      severe = weighted.mean(severe, person_days) * 365/12,
      time = mean(time),
      .by = c(age_lower, age_upper, year, month)
    ) %>%
    group_by(age_lower, age_upper) %>%
    mutate(cumul_clinical = cumsum(clinical)) %>% ungroup()

  infectivity_monthly <- sim %>%
    postie:::format_time(baseline_year = 2000) %>%
    select(starts_with('infectiv'), year, month) %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), mean),
      .by = c('month', 'year')
    )

  out <- left_join(inc, proc, by = c('year', 'month', 'time')) %>%
    left_join(infectivity_monthly, by = c('year','month')) %>%
    # add col for timestep of vaccination
    mutate(tbv_timing = if (tbv_coverage > 0) tbv_timestep_offset - burnin else NA_real_,
           draw = draw)

  return(out)
}

# Running across a range of parameter sets----

# Create a data frame of parameters
draws <- sample(1000, 20) # get 10 random draws from parameter sets
seasonalities <- c('seasonal','perennial')
tbv_coverages <- c(0.9)
tbv_ages_lab <- c('5+', 'SAC')
models <- c('m0', 'm1', 'm2')
pop <- 20000
combo <- crossing(seasonalities, tbv_coverages, tbv_ages_lab, models, pop, draws)
baseline <- data.frame(seasonalities = c('perennial','seasonal'),
                       tbv_coverages = 0,
                       tbv_ages_lab = 'noTBV',
                       models = 'baseline',
                       pop = pop) %>%
  crossing(draws)
combos <- rbind(combo, baseline)
combo_lst <- split(combos, seq_len(nrow(combos)))

# Map age labels to the ages passed to set_tbv (SAC = 5-15)
age_lookup <- list('5+' = 5:200, 'SAC' = 5:15)

cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
on.exit(parallel::stopCluster(cl))

parallel::clusterEvalQ(cl, {
  library(malariasimulation); library(postie); library(tidyverse)
})
parallel::clusterExport(cl, c('run_wtbaparams', 'model_pars', 'age_lookup',
                              'year'))  # plus anything else the function uses

results <- parallel::parLapply(cl, combo_lst, function(x){
  pars <- model_pars[model_pars$model == x$models, ]
  ages <- if (x$tbv_ages_lab == 'noTBV') 1 else age_lookup[[x$tbv_ages_lab]]

  out <- run_wtbaparams(
    m = pars$m,
    k = pars$k,
    pop = x$pop,
    draw = x$draws,
    season = x$seasonalities,
    tbv_coverage = x$tbv_coverages,
    tbv_ages = age_lookup[[x$tbv_ages_lab]]
  )

  message('finished ', x$seasonalities, ' | ', x$tbv_ages_lab, ' | ',
          x$models, ' | draw ', x$draw)

  out %>%
    mutate(seasonality = x$seasonalities,
           tbv_coverage = x$tbv_coverages,
           tbv_ages = x$tbv_ages_lab,
           model = x$models)
  })


all_res <- bind_rows(results) %>%
  filter(year < 2005 & year > 2000) %>%
  mutate(scenario = paste(seasonality, paste0('cov', tbv_coverage), tbv_ages, sep = '_'))

# 2. Cases averted vs m0 within each scenario ---------------------------------
join_keys <- c('draw', 'scenario', 'age_lower', 'age_upper', 'year', 'month')

m0_ref <- all_res %>%
  filter(model == 'm0') %>%
  select(all_of(join_keys), clinical_m0 = clinical, severe_m0 = severe)

tbv_res <- all_res %>% filter(model != 'baseline')

full <- all_res %>%
  left_join(m0_ref, by = join_keys) %>%
  mutate(ca_clinical = clinical_m0 - clinical,
         ca_severe = severe_m0 - severe) %>%
  arrange(scenario, model, draw, age_lower, time) %>%
  group_by(draw, scenario, model, age_lower) %>%
  mutate(cumul_ca_clinical = cumsum(ca_clinical)) %>%
  ungroup() %>%
  select(draw, scenario, seasonality, tbv_coverage, tbv_ages, model, tbv_timing,
         age_lower, age_upper, year, month, time,
         lm_prevalence_2_10, infectivity,
         clinical, clinical_m0, ca_clinical, cumul_ca_clinical,
         severe, severe_m0, ca_severe)

# 95% CrI across draws
summarise_cri <- function(df, vars, by){
  df %>%
    group_by(across(all_of(by))) %>%
    summarise(across(all_of(vars),
                     list(med = ~ median(.x, na.rm = TRUE),
                          lower  = ~ quantile(.x, 0.025, na.rm = TRUE, names = FALSE),
                          upper  = ~ quantile(.x, 0.975, na.rm = TRUE, names = FALSE)),
                     .names = '{.col}_{.fn}'),
              .groups = 'drop')
}


vars <- c('clinical', 'severe', 'lm_prevalence_2_10', 'infectivity',
          'ca_clinical', 'ca_severe', 'cumul_ca_clinical')
sum_all <- summarise_cri(full, vars,
                         by = c('scenario', 'seasonality', 'tbv_ages', 'model', 'tbv_timing',
                                'age_lower', 'age_upper', 'year', 'month', 'time'))
timing <- full %>% distinct(scenario, tbv_timing)

# 3. Plots for each scenario ----------------------------------------------------
model_cols_named <- setNames(model_cols, c('m0', 'm1', 'm2'))

make_plots <- function(df,
                       title){ # tbv timestep is this if the vaccine was delivered

  scales <- list(scale_color_manual(values = model_cols_named),
                 scale_fill_manual(values = model_cols_named))

  a_leg <- ggplot(df) +
    geom_line(aes(x = time, y = lm_prevalence_2_10_med, color = model)) +
    geom_ribbon(aes(x = time, ymin = lm_prevalence_2_10_lower, ymax = lm_prevalence_2_10_upper, fill = model),
                alpha = 0.2) +
    geom_vline(aes(xintercept = tbv_timing / 365 + 2000), color = 'red', linetype = 2) +
    scales +
    labs(x = 'Year', y = 'LM PfPR2-10') +
    kmtheme + theme(legend.position = 'bottom')
  a <- a_leg + theme(legend.position = 'none')

  b <- ggplot(df) +
    geom_line(aes(x = time, y = clinical_med, color = model)) +
    geom_ribbon(aes(x = time, ymin = clinical_lower, ymax = clinical_upper, fill = model),
                alpha = 0.2) +
    geom_vline(aes(xintercept = tbv_timing / 365 + 2000), color = 'red', linetype = 2) +
    facet_wrap(~age_lower, nrow = 3, scales = 'free_y') +
    scales +
    labs(x = 'Year', y = 'Monthly clinical incidence') +
    kmtheme + theme(legend.position = 'none')

  d <- ggplot(df) +
    geom_line(aes(x = time, y = ca_clinical_med, color = model)) +
    geom_ribbon(aes(x = time, ymin = ca_clinical_lower, ymax = ca_clinical_upper, fill = model),
                alpha = 0.2) +
    geom_vline(aes(xintercept = tbv_timing / 365 + 2000), color = 'red', linetype = 2) +
    geom_hline(aes(yintercept = 0), color = 'darkgreen', linetype = 2) +
    facet_wrap(~age_lower, nrow = 3, scales = 'free_y') +
    scales +
    labs(x = 'Year', y = 'Monthly clinical incidence averted (vs m0)') +
    kmtheme + theme(legend.position = 'none')

  e <- ggplot(df) +
    geom_line(aes(x = time, y = infectivity_med, color = model)) +
    geom_ribbon(aes(x = time, ymin = infectivity_lower, ymax = infectivity_upper, fill = model),
                alpha = 0.2) +
    geom_vline(aes(xintercept = tbv_timing / 365 + 2000), color = 'red', linetype = 2) +
    scales +
    labs(x = 'Year', y = 'Infectivity') +
    kmtheme + theme(legend.position = 'none')

  legenda <- cowplot::get_legend(a_leg)
  pp <- (b | d | (a / e)) + patchwork::plot_annotation(title = title)
  cowplot::plot_grid(pp, legenda, nrow = 2, rel_heights = c(1, 0.1))
}

out_dir <- 'stan_fitting/outputs/scenarios'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- unique(sum_all$scenario)
scenarios <- scenarios[-grep('noTBV', scenarios)]
plots <- lapply(setNames(scenarios, scenarios), function(s){
  p <- make_plots(df = filter(sum_all, scenario == s),
                  title = s)
  ggsave(file.path(out_dir, paste0('malariasim_n20k_', s, '20draws.png')),
         p, width = 12, height = 8)
  p
})
plots$seasonal_cov0.9_SAC
plots$`seasonal_cov0.9_5+`
plots$perennial_cov0.9_SAC
plots$`perennial_cov0.9_5+`

saveRDS(all_res, 'stan_fitting/outputs/scenarios/all_results_pop20k_202061019_20draws.rds')


# Run examples  ----
proc_orig <- run_wtbaparams(m = mALL, k = rrALL, pop = 10000) %>%
  mutate(model = 'm0')
proc_m1 <- run_wtbaparams(m = 0.3587039, k = 0.02142727, pop = 10000) %>% # from tbaplus_fitting.R
  mutate(model = 'm1')
proc_m2 <- run_wtbaparams(m = 0.0001570028, k = 2.049708e-05, pop = 10000) %>% # from tbaplus_fitting.R
  mutate(model = 'm2')

proc <- rbind(proc_orig, proc_m1, proc_m2) %>%
  filter(year < 2005 & year > 2000)

cases_averted <- proc %>% ungroup() %>%
  dplyr::select(age_lower, age_upper, clinical, model, time) %>%
  dplyr::mutate(dplyr::across(c(age_lower, age_upper),
                              ~ round(.x, digits = 0))) %>%
  pivot_wider(names_from = c('model'),
              values_from = 'clinical') %>%
  mutate(ca_m1 = m0 - m1,
         ca_m2 = m0 - m2)

# Plot prevalence 2-10 per person per year, averaged by month
a_leg <- ggplot(proc) +
  geom_line(aes(x = time, y = lm_prevalence_2_10, color = model)) +
  geom_vline(aes(xintercept = 2000+2), color = 'red', linetype = 2) +
  scale_color_manual(values = model_cols) +
  labs(x = 'Year',
       y = 'LM PfPR2-10') +
  kmtheme  +
  theme(legend.position = 'bottom')

a <- a_leg +
  theme(legend.position = 'none')


# Plot monthly incidence by age group
b <- ggplot(proc) +
  geom_line(aes(x = time, y = clinical, color = model)) +
  geom_vline(aes(xintercept = 2000+2), color = 'red', linetype = 2) +
  facet_wrap(~age_lower, nrow = 3,
             scales = 'free_y')+
  scale_color_manual(values = model_cols) +
  labs(x = 'Year',
       y = 'Monthly clinical incidence') +
  kmtheme +
  theme(legend.position = 'none')


# Plot monthly cumulative incidence by age group
c <- ggplot(proc) +
  geom_line(aes(x = time, y = cumul_clinical, color = model)) +
  geom_vline(aes(xintercept = 2000+2), color = 'red', linetype = 2) +
  facet_wrap(~age_lower, nrow = 3, scales = 'free')+
  scale_color_manual(values = model_cols) +
  kmtheme +
  theme(legend.position = 'none')


# Plot incidence averted by age group
d_leg <- ggplot(cases_averted %>% filter(year < 2005)) +
  geom_line(aes(x = time, y = ca_m1, color = 'm1')) +
  geom_line(aes(x = time, y = ca_m2, color = 'm2')) +
  geom_vline(aes(xintercept = 2000+2), color = 'red', linetype = 2) +
  geom_hline(aes(yintercept = 0), color = 'darkgreen', linetype = 2) +
  labs(y = 'Monthly clinical incidence averted (vs m0)',
       color = 'model',
       x = 'Year')+
  scale_color_manual(values = model_cols) +
  facet_wrap(~age_lower,
             nrow = 3,
             scales = 'free_y') +
  kmtheme+
  theme(legend.position = 'bottom')
d <- d_leg +
  theme(legend.position = 'none')

# Plot infectivity
e <- ggplot(proc) +
  geom_line(aes(x = time, y = infectivity, color = model))+
  geom_vline(aes(xintercept = 2000+2), color = 'red', linetype = 2) +
  scale_color_manual(values = model_cols) +
  labs(x = 'Year',
       y = 'Infectivity') +
  kmtheme +
  theme(legend.position = 'none')


library(patchwork)
legenda <- cowplot::get_legend(a_leg)
# legendd <- cowplot::get_legend(d_leg)
pp <- (b | d | (a / e))
# legends <- cowplot::plot_grid(legenda, legendd)
cowplot::plot_grid(pp, legenda, nrow = 2, rel_heights = c(1,0.1))
ggsave('stan_fitting/outputs/malariasim_n100k.png',
       width = 12, height = 8)



