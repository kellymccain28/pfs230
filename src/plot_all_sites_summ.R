# Script to plot infectivity for all sites, using summarised datasets over parameter draws
library(ggplot2)
library(tidyverse)

# Model outputs
# path <- 'outputs/2026-09-08_unweighted/' # first round of new calibration with Jen's data - no weighting
# path <- 'outputs/weighted_calibration/' # second round of new calibration with Jen's data - 1/SE
# path <- 'outputs/weighted_calibrationvariance/' # second round of new calibration with Jen's data - 1/SE^2
# path <- 'outputs/weighted_calibrationvariance_MAP/' # third round of new calibration with Jen's data - 1/SE^2 + site file MAP pfpr

make_plots <- function(path, label){
  all_summ <- readRDS(paste0(path, 'model_outputs/processed_output_summ.rds'))

  # Jen's data
  ento <- readRDS('data/Ento-reports-Aug2026/entomology_summarized.rds')
  pfpr <- readRDS('data/Ento-reports-Aug2026/parasitaemia_summarized.rds')


  # Get the monthly data and make long
  infectivity_monthly_summ <- all_summ$infectivity_monthly_summ %>%
    filter(year >= 2025) %>%
    select(country, site_name, ur, target_type, date, month, year,
           starts_with('mean_inf'), starts_with('prop_sum')
    ) %>%
    pivot_longer(
      cols = -c('country', 'site_name', 'ur', 'target_type', 'date', 'month', 'year'),
      names_to = c("metric", "age_group", "estimate"),
      names_pattern = "^(infectivity|prop_sum_inf|mean_inf)(?:_(under5|SAC|16plus|youngSAC|oldSAC|18plus))?_(median|lower|upper)$",
      values_to = 'value'
    ) %>%
    pivot_wider(
      names_from = c(metric, estimate),
      values_from = value,
      names_glue = "{metric}_{estimate}"
    ) %>%
    mutate(age_group = factor(age_group, levels = c('under5', 'SAC', '16plus',
                                                    'youngSAC','oldSAC', '18plus')))

  three_groups <- c('under5', 'SAC', '16plus')
  four_groups <- c('under5', 'youngSAC','oldSAC', '18plus')

  ltc_cols_type <- ltc::palettes$casa_natal
  ltc_cols_age <- ltc::palettes$expevo
  colors <- c('Data' = ltc_cols_type[1],
              'Model' = ltc_cols_type[3],
              'Simulated DSF positivity' = ltc_cols_type[2])


  # Plot the mean per person infectivity by age group
  p1 <- ggplot(infectivity_monthly_summ %>% filter(age_group %in% four_groups)) +
    geom_ribbon(aes(x = date, ymin = mean_inf_lower, ymax = mean_inf_upper,
                    group = age_group, fill = age_group),
                alpha = 0.3) +
    geom_line(aes(x = date, y = mean_inf_median, group = age_group, color = age_group)) +
    scale_color_manual(values = ltc_cols_age) +
    scale_fill_manual(values = ltc_cols_age) +
    scale_x_date(labels = scales::label_date_short()) +
    facet_wrap(~country, scales = 'free') +
    labs(x = 'Date',
         y = 'Per-bite (person) infectivity',
         color = 'Age group',
         fill = 'Age group') +
    theme_classic(base_size = 12)

  # Compare against data
  # For this, I want to simulate the DSFs using binomial draws with p = infectivity at each month in model output, matching ns in trial data
  all <- readRDS(paste0(path, 'all_processed_output.rds'))
  infect_monthly <- all %>%
    map('infectivity_monthly') %>%
    list_rbind() %>%
    mutate(date = floor_date(date, unit = 'month')) %>%
    filter(date %in% ento$date)

  df <- infect_monthly %>%
    select(date, country, site_name, parameter_draw,
           starts_with('mean_inf')) %>%
    pivot_longer(
      cols = -c(date, site_name, country, parameter_draw),
      names_to = c(".value", "age_group"),
      names_pattern = "(mean_inf)_(under5|SAC|16plus|youngSAC|oldSAC|18plus)"
    ) %>%
    filter(age_group =='SAC') %>%
    mutate(age_group = factor(age_group, levels = c('under5','SAC','16plus',
                                                    'youngSAC','oldSAC', '18plus'))) %>%
    left_join(ento %>%
                filter(age_group == '5-17') %>%
                mutate(age_group = 'SAC'), by = c('country','date', 'age_group')) %>%
    select(date, country, site_name, parameter_draw, age_group,
           n_positive_mosq, n_mosq_dissected, starts_with('mosq_positivity_rate'),
           starts_with('mean_inf'))

  df_1sim <- df %>%
    group_by(date, country, site_name, parameter_draw, age_group) %>%
    mutate(n_avg_dissected = mean(n_mosq_dissected)) %>%
    # Simulate 1 DSFs with ns from data and ps from model output
    rowwise() %>%
    mutate(sim_pos_dsf = rbinom(1, size = n_avg_dissected, prob = mean_inf),
           sim_mosq_pos_rate = sim_pos_dsf / n_avg_dissected) %>%
    ungroup()

  # ggplot(df_1sim) +
  #   geom_line(aes(x = date, y = mean_inf, group = parameter_draw, color = 'Model'), alpha = 0.5) +
  #   geom_line(aes(x = date, y = sim_mosq_pos_rate, group = parameter_draw, color = 'Simulated DSF positivity'), alpha = 0.7) +
  #   geom_pointrange(aes(x = date, y = mosq_positivity_rate/100,
  #                       ymin = mosq_positivity_rate_lower/100, ymax = mosq_positivity_rate_upper/100, color = 'Data'),
  #                   size = 0.3) +
  #   scale_color_manual(values = colors) +
  #   facet_wrap(~country, scales = 'free') +
  #   labs(x = 'Date',
  #        y = 'Per-person infectivity',
  #        color = NULL) +
  #   theme_classic(base_size = 12)

  # Simulate 100 DSFs per row, using n from data and p from model output
  df_100sims <- df %>%
    group_by(date, country, site_name, parameter_draw, age_group) %>%
    mutate(n_avg_dissected = mean(n_mosq_dissected)) %>%
    mutate(sim_pos_dsf = map2(n_avg_dissected, mean_inf,
                              ~ rbinom(100, size = .x, prob = .y))) %>%
    unnest(sim_pos_dsf)  %>%
    group_by(date, country, site_name, parameter_draw, age_group) %>%
    mutate(sim_id = row_number()) %>%
    ungroup() %>%
    mutate(sim_mosq_pos_rate = sim_pos_dsf / n_avg_dissected)

  # get summary over the simulated DSFs
  sim100_summary <- df_100sims %>%
    group_by(date, site_name, country) %>%
    summarise(across(c(sim_mosq_pos_rate, mean_inf),
                     .fns = list(
                       median = ~median(.x, na.rm = TRUE),
                       lower = ~quantile(.x, 0.025, na.rm = TRUE),
                       upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                     .names = "{.col}_{.fn}"),
              .groups = 'drop') %>%
    left_join(ento %>%
                select(date, age_group, country, mosq_positivity_rate, mosq_positivity_rate_lower, mosq_positivity_rate_upper) %>%
                filter(age_group == '5-17'))

  p2 <- ggplot(sim100_summary) +
    # simulated DSFs
    geom_line(aes(x = date, y = sim_mosq_pos_rate_median, color = 'Simulated DSF positivity'),
              alpha = 0.7, linewidth = 0.8) +
    geom_ribbon(aes(x = date, ymin = sim_mosq_pos_rate_lower, ymax = sim_mosq_pos_rate_upper,
                    fill = 'Simulated DSF positivity'), alpha = 0.3) +
    # model output
    geom_line(aes(x = date, y = mean_inf_median, color = 'Model'),
              alpha = 0.7, linewidth = 0.8) +
    geom_ribbon(aes(x = date, ymin = mean_inf_lower, ymax = mean_inf_upper,
                    fill = 'Model'), alpha = 0.3) +
    # # DSFs from data
    geom_pointrange(aes(x = date, y = mosq_positivity_rate/100,
                        ymin = mosq_positivity_rate_lower/100, ymax = mosq_positivity_rate_upper/100, color = 'Data'),
                    size = 0.3) +
    facet_wrap(~country, scales = 'free') +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_x_date(labels = scales::label_date_short()) +
    labs(x = 'Date',
         y = 'Per-person infectivity',
         color = NULL, fill = NULL) +
    theme_classic(base_size = 12)


  # Plot the relative infectivity by age group (proportion of summed population infectivity)
  p3a <- ggplot(infectivity_monthly_summ %>% filter(age_group %in% three_groups)) +
    geom_ribbon(aes(x = date, ymin = prop_sum_inf_lower, ymax = prop_sum_inf_upper,
                    group = age_group, fill = age_group),
                alpha = 0.3) +
    geom_line(aes(x = date, y = prop_sum_inf_median, group = age_group, color = age_group)) +
    scale_color_manual(values = ltc_cols_age) +
    scale_fill_manual(values = ltc_cols_age) +
    scale_x_date(labels = scales::label_date_short()) +
    facet_wrap(~country, scales = 'free') +
    labs(x = 'Date',
         y = 'Relative infectivity',
         color = 'Age group',
         fill = 'Age group') +
    theme_classic(base_size = 12)

  p3b <- ggplot(infectivity_monthly_summ %>% filter(age_group %in% four_groups)) +
    geom_ribbon(aes(x = date, ymin = prop_sum_inf_lower, ymax = prop_sum_inf_upper,
                    group = age_group, fill = age_group),
                alpha = 0.3) +
    geom_line(aes(x = date, y = prop_sum_inf_median, group = age_group, color = age_group)) +
    scale_color_manual(values = ltc_cols_age) +
    scale_fill_manual(values = ltc_cols_age) +
    scale_x_date(labels = scales::label_date_short()) +
    facet_wrap(~country, scales = 'free') +
    labs(x = 'Date',
         y = 'Relative infectivity',
         color = 'Age group',
         fill = 'Age group') +
    theme_classic(base_size = 12)

  # Summarize infectivity over last year of sim
  infectivity_lastyear_tbl <- all_summ$infectivity_annual_summ %>%
    filter(year == 2026) %>%
    select(country, site_name, ur, target_type, date, month, year,
           starts_with('mean_inf'), starts_with('prop_sum')
    ) %>%
    pivot_longer(
      cols = -c('country', 'site_name', 'ur', 'target_type', 'date', 'month', 'year'),
      names_to = c("metric", "age_group", "estimate"),
      names_pattern = "^(infectivity|prop_sum_inf|prop_mean_inf|mean_inf)(?:_(under5|SAC|16plus|youngSAC|oldSAC|18plus))?_(median|lower|upper)$",
      values_to = 'value'
    ) %>%
    pivot_wider(
      names_from = c(metric, estimate),
      values_from = value,
      names_glue = "{metric}_{estimate}"
    ) %>%
    mutate(age_group = factor(age_group, levels = c('under5', 'SAC', '16plus',
                                                    'youngSAC','oldSAC', '18plus')))

  p4a <- ggplot(infectivity_lastyear_tbl %>% filter(age_group %in% three_groups)) +
    geom_col(aes(x = age_group, y = prop_sum_inf_median, fill = age_group, color = age_group), alpha = 0.6) +
    geom_errorbar(aes(x = age_group, ymin = prop_sum_inf_lower, ymax = prop_sum_inf_upper, color = age_group),
                  width = 0.2, linewidth = 0.8) +
    geom_text(aes(x = age_group, y = 0.03, label = round(prop_sum_inf_median,2)),
              size = 3) +
    labs(y = 'Proportion of sum infectivity',
         x = 'Age group',
         fill = 'Age group',
         color = 'Age group')+
    scale_fill_manual(values = ltc_cols_age) +
    scale_color_manual(values = ltc_cols_age) +
    theme_classic(base_size = 12)  +
    facet_wrap(~country)

  p4b <- ggplot(infectivity_lastyear_tbl %>% filter(age_group %in% four_groups)) +
    geom_col(aes(x = age_group, y = prop_sum_inf_median, fill = age_group, color = age_group), alpha = 0.6) +
    geom_errorbar(aes(x = age_group, ymin = prop_sum_inf_lower, ymax = prop_sum_inf_upper, color = age_group),
                  width = 0.2, linewidth = 0.8) +
    geom_text(aes(x = age_group, y = 0.03, label = round(prop_sum_inf_median,2)),
              size = 3) +
    labs(y = 'Proportion of sum infectivity',
         x = 'Age group',
         fill = 'Age group',
         color = 'Age group')+
    scale_fill_manual(values = ltc_cols_age) +
    scale_color_manual(values = ltc_cols_age) +
    theme_classic(base_size = 12)  +
    facet_wrap(~country)

  p5a <- ggplot(infectivity_lastyear_tbl %>% filter(age_group %in% three_groups)) +
    geom_col(aes(x = country, y = prop_sum_inf_median, fill = age_group, color= age_group),
             position = position_fill(), alpha = 0.6)+
    geom_text(aes(x = country, y = prop_sum_inf_median, label = round(prop_sum_inf_median,2), group = age_group),
              position = position_fill(vjust = 0.5), size = 3) +
    labs(y = 'Proportion of sum infectivity',
         x = 'Trial site',
         fill = 'Age group',
         color = 'Age group') +
    scale_fill_manual(values = ltc_cols_age) +
    scale_color_manual(values = ltc_cols_age) +
    theme_classic(base_size = 12)

  p5b <- ggplot(infectivity_lastyear_tbl %>% filter(age_group %in% four_groups)) +
    geom_col(aes(x = country, y = prop_sum_inf_median, fill = age_group, color= age_group),
             position = position_fill(), alpha = 0.6)+
    geom_text(aes(x = country, y = prop_sum_inf_median+0.05, label = round(prop_sum_inf_median,2), group = age_group),
              position = position_fill(vjust = 0.5), size = 3) +
    labs(y = 'Proportion of sum infectivity',
         x = 'Trial site',
         fill = 'Age group',
         color = 'Age group') +
    scale_fill_manual(values = ltc_cols_age) +
    scale_color_manual(values = ltc_cols_age) +
    theme_classic(base_size = 12)

  p6a <- ggplot(infectivity_lastyear_tbl %>% filter(age_group %in% three_groups)) +
    geom_col(aes(x = age_group, y = mean_inf_median, fill = age_group, color = age_group), alpha = 0.6) +
    geom_errorbar(aes(x = age_group, ymin = mean_inf_lower, ymax = mean_inf_upper, color = age_group),
                  width = 0.2, linewidth = 0.8) +
    geom_text(aes(x = age_group, y = 0.001, label = round(prop_sum_inf_median,2)),
              size = 3) +
    labs(y = 'Mean per-person infectivity',
         x = 'Age group',
         fill = 'Age group',
         color = 'Age group')+
    scale_fill_manual(values = ltc_cols_age) +
    scale_color_manual(values = ltc_cols_age) +
    theme_classic(base_size = 12)  +
    facet_wrap(~country)

  p6b <- ggplot(infectivity_lastyear_tbl %>% filter(age_group %in% four_groups)) +
    geom_col(aes(x = age_group, y = mean_inf_median, fill = age_group, color = age_group), alpha = 0.6) +
    geom_errorbar(aes(x = age_group, ymin = mean_inf_lower, ymax = mean_inf_upper, color = age_group),
                  width = 0.2, linewidth = 0.8) +
    geom_text(aes(x = age_group, y = 0.001, label = round(prop_sum_inf_median,2)),
              size = 3) +
    labs(y = 'Mean per-person infectivity',
         x = 'Age group',
         fill = 'Age group',
         color = 'Age group')+
    scale_fill_manual(values = ltc_cols_age) +
    scale_color_manual(values = ltc_cols_age) +
    theme_classic(base_size = 12)  +
    facet_wrap(~country)

  #####################################################
  ### Plot epi outputs
  #####################################################
  monthly_epi <- all %>%
    map('monthly_epi_output') %>%
    list_rbind()

  # long datasets
  prev_monthly <- monthly_epi %>%
    dplyr::select(time, contains('prev'), month, ur, country, parameter_draw,
                  site_name, target_type) %>%
    pivot_longer(cols = contains('prev'),
                 names_to = c(".value", "age_group_prev"),
                 names_pattern = "(lm_prevalence)_(.*)") %>%
    distinct() %>%
    mutate(date = as.Date(format(floor_date(date_decimal(time)), '%Y-%m-%d')),
           date = floor_date(date, unit ='month')) %>%
    filter(age_group_prev %in% c('5_9','9_18')) %>%
    mutate(age_group = case_when(
      age_group_prev == '5_9' ~ '5-8',
      age_group_prev == '9_18' ~ '9-17',
      TRUE ~ NA
    )) %>%
    filter(time >= 2025) %>%
    group_by(time, age_group, month, ur, country,
             site_name, target_type) %>%
    mutate(lm_prevalence_med = median(lm_prevalence))

  p7 <- ggplot(prev_monthly) +
    geom_line(aes(x = date, y = lm_prevalence, group = parameter_draw, color = age_group), alpha= 0.1) +
    geom_line(aes(x = date, y = lm_prevalence_med, color = age_group), linewidth = 0.5, alpha = 0.8) +
    geom_pointrange(data = pfpr %>% filter(age_group == '9-17' | age_group == '5-8'),
                    aes(x = date,
                        y = pf_positivity_rate/100,
                        ymin = pf_positivity_rate_lower/100,
                        ymax = pf_positivity_rate_upper/100,
                        group = age_group,
                        color = age_group),
                    position = position_dodge(width = 20),
                    size = 0.3) +
    facet_grid(country ~ age_group, scales = 'free')+
    scale_color_manual(values = ltc_cols_type) +
    scale_x_date(labels = scales::label_date_short()) +
    labs(x = 'Time',
         y = 'LM PfPR',
         color = NULL) +
    theme_classic(base_size = 12)

  ggsave(paste0(path, label, "_infectivity_mean_perperson.pdf"), p1, height = 5, width = 10)
  ggsave(paste0(path, label, "_infectivity_mean_perperson_data_100sims.pdf"), p2, height = 5, width = 10)
  ggsave(paste0(path, label, "_relative_infectivity_byage_overtime_3ages.pdf"), p3a, height = 5, width = 10)
  ggsave(paste0(path, label, "_relative_infectivity_byage_overtime_4ages.pdf"), p3b, height = 5, width = 10)
  ggsave(paste0(path, label, "_relative_infectivity_byage_lastyear_3ages.pdf"), p4a)
  ggsave(paste0(path, label, "_relative_infectivity_byage_lastyear_4ages.pdf"), p4b)
  ggsave(paste0(path, label, "_relative_infectivity_byage_lastyear_posfill_3ages.pdf"), p5a)
  ggsave(paste0(path, label, "_relative_infectivity_byage_lastyear_posfill_4ages.pdf"), p5b)
  ggsave(paste0(path, label, "_infectivity_mean_perperson_byage_lastyear_3ages.pdf"), p6a)
  ggsave(paste0(path, label, "_infectivity_mean_perperson_byage_lastyear_4ages.pdf"), p6b)
  ggsave(paste0(path, label, "_prevalence_modeltodata.pdf"), p7, height = 8, width = 10)

}

# BElow won't work anymore because they don't have new infectivity age groups
# make_plots(path = 'outputs/2026-09-08_unweighted/', label = 'unweighted')# first round of new calibration with Jen's data - no weighting
# make_plots(path = 'outputs/weighted_calibration/', label = 'weighted_SE')# second round of new calibration with Jen's data - 1/SE
# make_plots(path = 'outputs/weighted_calibrationvariance/', label = 'weighted_var')# second round of new calibration with Jen's data - 1/SE^2
# make_plots(path = 'outputs/weighted_calibrationvariance_MAP/', label = 'weighted_var_MAP')# third round of new calibration with Jen's data - 1/SE^2 + MAP pfpr

make_plots(path = 'outputs/jen_map_test_new_infectivity_cluster/', label = 'testing')
