# Script to run different scenarios assuming different models for the relationship between oocysts --> P(inf)
# which lead to different relationships between TRA-TBA+
# pak::pak("mrc-ide/malariasimulation")
library(malariasimulation) # make sure this is the version on the main branch
library(postie)
library(ggplot2)
library(tidyverse)

model_cols <- c('#ff595e','#8ac926','#6a4c93')
kmtheme <- theme_bw(base_size = 12)
# Bompard parameters for oocyst distribtuion (m0)
mALL <- 0.000157   # NB mean
rrALL <- 0.00000495 # NB dispersion (size)

year = 365

# original model parameterss
m0_pars <- c(model = 'm0', m = mALL, k = rrALL)
m1_pars <- c(model = 'm1', m = 0.3587039, k = 0.02142727)
m2_pars <- c(model = 'm2', m = 0.0001570028, k = 2.049708e-05)
#updated parameters (4/11/2026) with rnbinom draws in loop
m0_pars <- c(model = 'm0', m = mALL, k = rrALL)
m1_pars <- c(model = 'm1', m = 6.902878e-06, k = 4.633258e-07)
m2_pars <- c(model = 'm2', m = 7.209415e-05, k = 1.077945e-05)
#updated parameters (9/11/2026) with solution dnbinom
m0_pars <- c(model = 'm0', m = mALL, k = rrALL)
m1_pars <- c(model = 'm1', m = 3.088107e-05, k = 1.648756e-06)
m2_pars <- c(model = 'm2', m = 8.257316e-05, k = 9.411515e-06)

model_pars <- data.frame(rbind(m0_pars, m1_pars, m2_pars)) %>%
  mutate(m = as.numeric(m),
         k = as.numeric(k))
rownames(model_pars) <- NULL

run_wtbaparams <- function(m, k,
                           pop,
                           draw,
                           season = 'highly_seasonal',
                           tbv_timestep = 12*year,
                           tbv_coverage = 0.8,
                           tbv_ages = 5:200,
                           timesteps = 20*year,
                           burnin = 10*year
                           ){

  if(timesteps < tbv_timestep){
    stop('TBV should be delivered before max timestep')}
  if(burnin > tbv_timestep){
    stop('TBV should be delivered after burnin')}

  # seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
  # drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
  # g0 = a0, a = g, b = h
  seas_name <- 'highly_seasonal'
  seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
  s1 <- tibble(seasonality, seas_name)

  seas_name <- 'seasonal'
  seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
  s2 <- tibble(seasonality, seas_name)

  seas_name <- 'perennial'
  seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
  s3 <- tibble(seasonality, seas_name)
  seasonality <- bind_rows(s1, s2, s3) %>%
    filter(seas_name == season)
  if (nrow(seasonality) != 1) stop('Unknown seasonality: ', season)

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
    set_clinical_treatment(1, 0.5, 0.5) #%>%
    # set_epi_outputs(
    #   age_group = c(0, 5, 15, 200) * year,
    #   prevalence = c(2, 10) * year,
    #   clinical_incidence = c(0, 5, 15, 200) * year,
    #   severe_incidence = c(0, 5, 15, 200) * year
    # )
  params$clinical_incidence_rendering_min_ages <- c(0, 5, 15, 0) * year
  params$clinical_incidence_rendering_max_ages <- c(5, 15, 200, 200) * year - 1
  params$severe_incidence_rendering_min_ages <- c(0, 5, 15, 0) * year
  params$severe_incidence_rendering_max_ages <- c(5, 15, 200, 200) * year - 1
  params$age_group_rendering_min_ages <- c(0, 5, 15, 0) * year
  params$age_group_rendering_max_ages <- c(5, 15, 200, 200) * year - 1
  params$prevalence_rendering_min_ages <- 2 * year
  params$prevalence_rendering_max_ages <- 10 * year - 1

  if(draw > 0) {
    params <- params %>%
      set_parameter_draw(draw)
  }

  # Set bednets
  # parameters below from Ellie's paper https://doi.org/10.1016/S2542-5196(21)00296-5
  # The parameter values for pyrethroid-only and pyrethroid-PBO nets at various resistance levels can be found in Table S1.3 in the Supplementary Appendix 2 of Sherrard-Smith, et al., 2022.
  r <- 0.609 # probability that mosquito repelled/blocked from feeding
  d <- 0.322 # probability that mosquito can be killed & prevented from feeding
  H <- 2.022 # half life of decay of a ITN's lethality
  meandur <- H/ log(2) # mean duration of effect. We convert the half life by dividing its value by log(2)

  # set timesteps every 3 years from 9 years before simulation
  itn_timesteps <- seq(burnin - 9*year, burnin + timesteps, 3*year)
  params <- set_bednets(params,
                        timesteps = itn_timesteps,
                        coverages = rep(0.6, length(itn_timesteps)),
                        retention = 4 * year,
                        dn0 = matrix(rep(d, length(itn_timesteps)), nrow = length(itn_timesteps), ncol = 1),
                        rn = matrix(rep(r, length(itn_timesteps)), nrow = length(itn_timesteps), ncol = 1),
                        rnm = matrix(rep(0.24, length(itn_timesteps)), nrow = length(itn_timesteps), ncol = 1),
                        gamman = rep(meandur, length(itn_timesteps)))

  params$phi_bednets <- 0.79 # to match Ellie's PNAS paper https://doi.org/10.1073/pnas.1820646116
  # "On average 79% of bites by the major malaria vectors occur during the time when people are in bed. This estimate is substantially lower than previous predictions, with results suggesting a nearly 10% lower proportion of bites taken at the time when people are beneath LLINs since the year 2000."

  # Parameterize the TBV
  peak <- peak_season_offset(params)
  tbv_timestep_offset <- tbv_timestep + peak - (3 * round(365/12))

  # primary series and booster dose
  params <- set_tbv(params,
                    timesteps = c(tbv_timestep_offset, tbv_timestep_offset + 365), # 3 months before peak and 1 year later
                    coverages = rep(tbv_coverage, 2),
                    ages = tbv_ages)

  # Set correlation between primary series and booster
  correlations <- get_correlation_parameters(params)
  correlations$inter_round_rho('tbv', 1)

  params$tbv_mu <- m
  params$tbv_mt <- m
  params$tbv_md <- m
  params$tbv_ma <- m
  params$tbv_k <- k

  # Set equilibrium
  params <- set_equilibrium(params, init_EIR = 10)

  # Run simulation with seed
  set.seed(draw) # to allow for paired comparisons between scenarios
  sim <- run_simulation(timesteps = timesteps,
                        parameters = params,
                        correlations = correlations)
  sim <- drop_burnin(sim, burnin = burnin)

  # Process raw output
  proc <- get_prevalence(sim) %>%
    summarise(
      lm_prevalence_2_10 = mean(lm_prevalence_2_10),
      time = mean(time),
      .by = c(year, month)
    )

  inc <- get_rates(sim) %>%
    arrange(age_lower, age_upper, time) %>%
    group_by(age_lower, age_upper) %>%
    mutate(cumul_clinical = cumsum(clinical) * 1000,   # per 1000 people
           cumul_severe   = cumsum(severe) * 1000) %>%
    ungroup()

  inc <- inc %>%
    summarise(
      cases_clinical = sum(clinical * person_days),,
      clinical = weighted.mean(clinical, person_days) * 365/12 * 1000, # per 1000 people
      severe = weighted.mean(severe, person_days) * 365/12 * 1000, # per 1000 people
      cumul_clinical = max(cumul_clinical),
      cumul_severe   = max(cumul_severe),
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

  doses <- sim %>%
    postie:::format_time(baseline_year = 2000) %>%
    dplyr::summarise(
      n_tbv_doses = sum(n_vaccinated_tbv, na.rm = TRUE),
      .by = c('month','year')
    )

  out <- left_join(inc, proc, by = c('year', 'month', 'time')) %>%
    left_join(infectivity_monthly, by = c('year','month')) %>%
    left_join(doses, by = c('year','month')) %>%
    # add col for timestep of vaccination
    mutate(tbv_timing = if (tbv_coverage > 0) tbv_timestep_offset - burnin else NA_real_,
           draw = draw,
           bednets = params$bednets)

  return(out)
}

# Running across a range of parameter sets----

# Create a data frame of parameters
draws <- 0 #sample(1000, 50) # get 10 random draws from parameter sets
seasonalities <- c('highly_seasonal')#'seasonal','perennial',
tbv_coverages <- c(0.8)
tbv_ages_lab <- c('5+', 'SAC')
models <- c('m0', 'm1', 'm2')
pop <- 100000
combo <- crossing(seasonalities, tbv_coverages, tbv_ages_lab, models, pop, draws)
baseline <- data.frame(seasonalities = 'highly_seasonal',#c('perennial','seasonal'),
                       tbv_coverages = 0,
                       tbv_ages_lab = 'noTBV',
                       models = 'm0',
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
          x$models, ' | draw ', x$draws)

  out %>%
    mutate(seasonality = x$seasonalities,
           tbv_coverage = x$tbv_coverages,
           tbv_ages = x$tbv_ages_lab,
           model = if (x$tbv_ages_lab == 'noTBV') 'baseline' else x$models)
  })


all_res <- bind_rows(results) %>%
  filter(year < 2005 & year > 2000) %>%
  mutate(scenario = paste(seasonality, paste0('cov', tbv_coverage), tbv_ages, sep = '_'))
saveRDS(all_res, 'stan_fitting/outputs/scenarios/all_results_pop100k_202061106_0draws_3mtbv_booster_newmodels.rds')

# 2. Cases averted vs m0 within each scenario ---------------------------------
join_keys <- c('draw', 'scenario', 'age_lower', 'age_upper', 'year', 'month')
baseline_keys <- c('draw', 'seasonality', 'age_lower', 'age_upper', 'year', 'month')

m0_ref <- all_res %>%
  filter(model == 'm0', tbv_ages != 'noTBV') %>%
  select(all_of(join_keys),
         cases_m0 = cases_clinical,
         clinical_m0 = clinical,
         severe_m0 = severe,
         cumul_clinical_m0 = cumul_clinical)

notbv_ref <- all_res %>%
  filter(tbv_ages == 'noTBV') %>% #same as model == 'baseline'
  select(all_of(baseline_keys),
         cases_baseline = cases_clinical,
         clinical_baseline = clinical,
         severe_baseline = severe,
         cumul_clinical_baseline = cumul_clinical)

full <- all_res %>% ungroup() %>%
  filter(model != 'baseline') %>%
  left_join(m0_ref, by = join_keys) %>%
  left_join(notbv_ref, by = baseline_keys) %>%
  arrange(scenario, model, draw, age_lower, age_upper, time) %>%
  group_by(draw, scenario, model, age_lower, age_upper) %>%
  mutate(ca_clinical_vsm0 = clinical_m0 - clinical,
         ca_clinical = clinical_baseline - clinical,
         cumul_ca_vsm0_clinical = cumul_clinical_m0 - cumul_clinical,
         cumul_ca_clinical = cumul_clinical_baseline - cumul_clinical,

         totaldoses = sum(n_tbv_doses),
         ca_cases            = cases_baseline - cases_clinical,
         ca_clinical_perdose = ca_cases / totaldoses * 1000,
         cumul_ca_clinical_perdose = cumsum(ca_cases) / totaldoses * 1000,

         ca_clinical_vsm0_perdose = ca_clinical_vsm0 / totaldoses * 1000,
         cumul_ca_vsm0_clinical_perdose = cumul_ca_vsm0_clinical / totaldoses * 1000,

         prop_ca_cumul = cumul_ca_clinical / cumul_clinical_baseline,  # cumulative proportion averted
         prop_ca_month = ca_clinical / clinical_baseline) %>%         # monthly proportion averted) %>%
  ungroup() %>%
  select(draw, scenario, seasonality, tbv_coverage, tbv_ages, model, tbv_timing,
         age_lower, age_upper, year, month, time,
         n_tbv_doses, totaldoses,
         ca_clinical_perdose, cumul_ca_clinical_perdose,
         ca_clinical_vsm0_perdose, cumul_ca_vsm0_clinical_perdose,
         lm_prevalence_2_10, infectivity,
         cases_clinical, cases_baseline,cases_m0,
         clinical, clinical_m0, clinical_baseline,
         ca_clinical_vsm0, ca_clinical, ca_cases,
         cumul_ca_clinical, cumul_ca_vsm0_clinical,
         cumul_clinical, cumul_clinical_m0, cumul_clinical_baseline,
         prop_ca_cumul, prop_ca_month)

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

vars <- c('clinical', 'clinical_baseline', 'clinical_m0',
          'n_tbv_doses','totaldoses',
          'ca_clinical_perdose', 'cumul_ca_clinical_perdose',
          'ca_clinical_vsm0_perdose', 'cumul_ca_vsm0_clinical_perdose',
          'lm_prevalence_2_10', 'infectivity',
          'cases_clinical','cases_baseline','cases_m0',
          'ca_cases', 'ca_clinical', 'ca_clinical_vsm0',
          'cumul_ca_clinical', 'cumul_ca_vsm0_clinical',
          'cumul_clinical', 'cumul_clinical_baseline', 'cumul_clinical_m0',
          'prop_ca_cumul', 'prop_ca_month')
sum_all <- summarise_cri(full, vars,
                         by = c('scenario', 'seasonality', 'tbv_ages', 'model', 'tbv_timing',
                                'age_lower', 'age_upper', 'year', 'month', 'time')) %>%
  mutate(age_upper = round(age_upper, 0),
         age_group = paste0(age_lower, '-', age_upper),
         age_group = factor(age_group, levels = c('0-5','5-15','15-200','0-200')))
timing <- full %>% distinct(scenario, tbv_timing)

saveRDS(sum_all, 'stan_fitting/outputs/scenarios/summarized_allresults_pop100k_20261106_0draws_3mtbv_booster_newmodels.rds')

# 3. Plots for each scenario ----------------------------------------------------
model_cols_named <- setNames(model_cols, c('m0', 'm1', 'm2'))

make_plots <- function(df,
                       title){ # tbv timestep is this if the vaccine was delivered

  scalesmod <- list(scale_color_manual(values = model_cols_named),
                 scale_fill_manual(values = model_cols_named))
  tbv_line <- geom_vline(aes(xintercept = tbv_timing / 365 + 2000),
                         color = 'red', linetype = 2)
  zero_line <- geom_hline(aes(yintercept = 0), color = 'darkgreen', linetype = 2)
  facet <- facet_wrap(~age_group, nrow = 3)

  # ---- shared panels ----
  a_leg <- ggplot(df) +
    geom_line(aes(x = time, y = lm_prevalence_2_10_med, color = model)) +
    # geom_ribbon(aes(x = time, ymin = lm_prevalence_2_10_lower, ymax = lm_prevalence_2_10_upper, fill = model),
    #             alpha = 0.2) +
    tbv_line +
    scalesmod +
    labs(x = 'Year', y = 'LM PfPR2-10') +
    kmtheme + theme(legend.position = 'bottom')
  a <- a_leg + theme(legend.position = 'none')

  e <- ggplot(df) +
    geom_line(aes(x = time, y = infectivity_med, color = model)) +
    # geom_ribbon(aes(x = time, ymin = infectivity_lower, ymax = infectivity_upper, fill = model),
    #             alpha = 0.2) +
    tbv_line +
    scalesmod +
    labs(x = 'Year', y = 'Infectivity') +
    kmtheme + theme(legend.position = 'none')

  # Monthly incidence, with the no-TBV counterfactual as a dashed black line
  b <- ggplot(df) +
    geom_line(aes(x = time, y = clinical_baseline_med),
              color = 'black', linetype = 3) +
    geom_line(aes(x = time, y = clinical_med, color = model)) +
    tbv_line +
    facet +# _wrap(~age_group, nrow = 3, scales = 'free_y') +
    scalesmod +
    labs(x = 'Year', y = 'Monthly clinical incidence\n(dotted = no TBV)') +
    kmtheme + theme(legend.position = 'none')

  # ---- Row 1: vs m0 ----
  d_m0 <- ggplot(df) +
    geom_line(aes(x = time, y = ca_clinical_vsm0_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    labs(x = 'Year', y = 'Monthly cases averted per 1000 people (vs m0)') +
    kmtheme + theme(legend.position = 'none')

  c_m0 <- ggplot(df) +
    geom_line(aes(x = time, y = cumul_ca_vsm0_clinical_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    labs(x = 'Year', y = 'Cumulative cases averted per 1000 people (vs m0)') +
    kmtheme + theme(legend.position = 'none')

  d_perdosem0 <- ggplot(df) +
    geom_line(aes(x = time, y = ca_clinical_vsm0_perdose_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    labs(x = 'Year', y = 'Monthly cases averted per person per 1000 doses (vs m0)') +
    kmtheme + theme(legend.position = 'none')

  c_perdosem0 <- ggplot(df) +
    geom_line(aes(x = time, y = cumul_ca_vsm0_clinical_perdose_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    labs(x = 'Year', y = 'Cumulative cases averted per person per 1000 doses (vs m0)') +
    kmtheme + theme(legend.position = 'none')

  # ---- Row 2: vs no TBV ----
  d_tbv <- ggplot(df) +
    geom_line(aes(x = time, y = ca_clinical_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    labs(x = 'Year', y = 'Monthly cases averted per 1000 people (vs no TBV)') +
    kmtheme + theme(legend.position = 'none')

  c_tbv <- ggplot(df) +
    geom_line(aes(x = time, y = cumul_ca_clinical_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    labs(x = 'Year', y = 'Cumulative cases averted per 1000 people (vs no TBV)') +
    kmtheme + theme(legend.position = 'none')

  d_perdosetbv <- ggplot(df) +
    geom_line(aes(x = time, y = ca_clinical_perdose_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    labs(x = 'Year', y = 'Monthly cases averted per person per 1000 doses (vs no TBV)') +
    kmtheme + theme(legend.position = 'none')

  c_perdosetbv <- ggplot(df) +
    geom_line(aes(x = time, y = cumul_ca_clinical_perdose_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    labs(x = 'Year', y = 'Cumulative cases averted per person per 1000 doses (vs no TBV)') +
    kmtheme + theme(legend.position = 'none')

  p_tbv <- ggplot(df) +
    geom_line(aes(x = time, y = prop_ca_cumul_med, color = model)) +
    tbv_line + zero_line +
    facet +
    scalesmod +
    scale_y_continuous(labels = scales::percent) +
    labs(x = 'Year', y = 'Cumulative proportion averted (vs no TBV)') +
    kmtheme + theme(legend.position = 'none')

  legenda <- cowplot::get_legend(a_leg)

  row1 <- (b | d_tbv | c_tbv | (a / e)) #d_m0 | c_m0 | d_perdosem0 | c_perdosem0 |
  row2 <- (d_perdosetbv | c_perdosetbv | p_tbv)

  pp <- (row1 / row2) + patchwork::plot_annotation(title = title)
  cowplot::plot_grid(pp, legenda, nrow = 2, rel_heights = c(1, 0.05))
}

out_dir <- 'stan_fitting/outputs/scenarios'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- unique(sum_all$scenario)
scenarios <- scenarios[!grepl('noTBV', scenarios)]
plots <- lapply(setNames(scenarios, scenarios), function(s){
  p <- make_plots(df = filter(sum_all, scenario == s,
                              (age_upper == 5 | age_upper == 200) & age_lower == 0),
                  title = s)
  ggsave(file.path(out_dir, paste0('malariasim_n100k_', s, '202061106_0draws_booster_newmodels.png')),
         p, width = 16, height = 16)
  p
})
plots$`highly_seasonal_cov0.8_5+`
plots$highly_seasonal_cov0.8_SAC
# plots$seasonal_cov0.9_SAC
# plots$`seasonal_cov0.9_5+`
# plots$perennial_cov0.9_SAC
# plots$`perennial_cov0.9_5+`



# Run examples  ----
proc_orig <- run_wtbaparams(m = mALL, k = rrALL, pop = 10000) %>%
  mutate(model = 'm0')
proc_m1 <- run_wtbaparams(m = 0.3587039, k = 0.02142727, pop = 10000) %>% # from tbaplus_fitting.R (when doing rnbinom() draws)
  mutate(model = 'm1')
proc_m2 <- run_wtbaparams(m = 0.0001570028, k = 2.049708e-05, pop = 10000) %>% # from tbaplus_fitting.R (when doing rnbinom() draws)
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



