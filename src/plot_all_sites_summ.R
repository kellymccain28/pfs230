# Script to plot infectivity for all sites, using summarised datasets over parameter draws
library(ggplot2)
library(tidyverse)

# Model outputs
path <- 'outputs/2026-09-08_weighted/' # first round of new calibration with Jen's data
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
    names_pattern = "^(infectivity|prop_sum_inf|prop_mean_inf|mean_inf)(?:_(under5|SAC|16plus))?_(median|lower|upper)$",
    values_to = 'value'
  ) %>%
  pivot_wider(
    names_from = c(metric, estimate),
    values_from = value,
    names_glue = "{metric}_{estimate}"
  ) %>%
  mutate(age_group = factor(age_group, levels = c('under5', 'SAC', '16plus')))

ltc_cols_type <- ltc::palettes$casa_natal
ltc_cols_age <- ltc::palettes$expevo
colors <- c('Data' = ltc_cols_type[1],
            'Model' = ltc_cols_type[3],
            'Simulated DSF positivity' = ltc_cols_type[2])


# Plot the mean per person infectivity by age group
ggplot(infectivity_monthly_summ) +
  geom_ribbon(aes(x = date, ymin = mean_inf_lower, ymax = mean_inf_upper,
                  group = age_group, fill = age_group),
              alpha = 0.3) +
  geom_line(aes(x = date, y = mean_inf_median, group = age_group, color = age_group)) +
  scale_color_manual(values = ltc_cols_age) +
  scale_fill_manual(values = ltc_cols_age) +
  facet_wrap(~site_name, scales = 'free') +
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
    names_pattern = "(mean_inf)_(under5|SAC|16plus)"
  ) %>%
  filter(age_group =='SAC') %>%
  mutate(age_group = factor(age_group, levels = c('under5','SAC','16plus'))) %>%
  left_join(ento %>%
              filter(age_group == '5-17') %>%
              mutate(age_group = 'SAC'), by = c('country','date', 'age_group')) %>%
  select(date, country, site_name, parameter_draw, age_group,
         n_positive_mosq, n_mosq_dissected, starts_with('mosq_positivity_rate'),
         starts_with('mean_inf'))

df_1sim <- df %>%
  # Simulate 1 DSFs with ns from data and ps from model output
  rowwise() %>%
  mutate(sim_pos_dsf = rbinom(1, size = n_mosq_dissected, prob = mean_inf),
         sim_mosq_pos_rate = sim_pos_dsf / n_mosq_dissected) %>%
  ungroup()

ggplot(df_1sim) +
  geom_line(aes(x = date, y = mean_inf, group = parameter_draw, color = 'Model'), alpha = 0.5) +
  geom_line(aes(x = date, y = sim_mosq_pos_rate, group = parameter_draw, color = 'Simulated DSF positivity'), alpha = 0.7) +
  geom_pointrange(aes(x = date, y = mosq_positivity_rate/100,
                      ymin = mosq_positivity_rate_lower/100, ymax = mosq_positivity_rate_upper/100, color = 'Data'),
                  size = 0.3) +
  scale_color_manual(values = colors) +
  facet_wrap(~site_name, scales = 'free') +
  labs(x = 'Date',
       y = 'Per-person infectivity',
       color = NULL) +
  theme_classic(base_size = 12)

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

ggplot(sim100_summary) +
  # model output
  geom_line(aes(x = date, y = mean_inf_median, group = parameter_draw, color = 'Model'),
            alpha = 0.7, linewidth = 0.8) +
  geom_ribbon(aes(x = date, ymin = mean_inf_lower, ymax = mean_inf_upper,
                  fill = 'Model'), alpha = 0.3) +
  # simulated DSFs
  geom_line(aes(x = date, y = sim_mosq_pos_rate_median, group = parameter_draw, color = 'Simulated DSF positivity'),
            alpha = 0.7, linewidth = 0.8) +
  geom_ribbon(aes(x = date, ymin = sim_mosq_pos_rate_lower, ymax = sim_mosq_pos_rate_upper,
                  fill = 'Simulated DSF positivity'), alpha = 0.3) +
  # # DSFs from data
  geom_pointrange(aes(x = date, y = mosq_positivity_rate/100,
                      ymin = mosq_positivity_rate_lower/100, ymax = mosq_positivity_rate_upper/100, color = 'Data'),
                  size = 0.3) +
  facet_wrap(~site_name, scales = 'free') +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = 'Date',
       y = 'Per-person infectivity',
       color = NULL, fill = NULL) +
  theme_classic(base_size = 12)


# Plot the relative infectivity by age group (proportion of summed population infectivity)
ggplot(infectivity_monthly_summ) +
  geom_ribbon(aes(x = date, ymin = prop_sum_inf_lower, ymax = prop_sum_inf_upper,
                  group = age_group, fill = age_group),
              alpha = 0.3) +
  geom_line(aes(x = date, y = prop_sum_inf_median, group = age_group, color = age_group)) +
  scale_color_manual(values = ltc_cols_age) +
  scale_fill_manual(values = ltc_cols_age) +
  facet_wrap(~site_name, scales = 'free') +
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
    names_pattern = "^(infectivity|prop_sum_inf|prop_mean_inf|mean_inf)(?:_(under5|SAC|16plus))?_(median|lower|upper)$",
    values_to = 'value'
  ) %>%
  pivot_wider(
    names_from = c(metric, estimate),
    values_from = value,
    names_glue = "{metric}_{estimate}"
  ) %>%
  mutate(age_group = factor(age_group, levels = c('under5', 'SAC', '16plus')))

ggplot(infectivity_lastyear_tbl) +
  geom_col(aes(x = age_group, y = prop_sum_inf_median, fill = age_group), alpha = 0.7) +
  geom_errorbar(aes(x = age_group, ymin = prop_sum_inf_lower, ymax = prop_sum_inf_upper, color = age_group),
                width = 0.2) +
  geom_text(aes(x = age_group, y = 0.03, label = round(prop_sum_inf_median,2)),
            size = 3) +
  labs(y = 'Proportion of sum infectivity',
       x = 'Age group',
       fill = 'Age group',
       color = 'Age group')+
  scale_fill_manual(values = ltc_cols_age) +
  scale_color_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12)  +
  facet_wrap(~site_name)

ggplot(infectivity_lastyear_tbl) +
  geom_col(aes(x = site_name, y = prop_sum_inf_median, fill = age_group),
           position = position_fill())+
  labs(y = 'Proportion of sum infectivity',
       x = 'Trial site',
       fill = 'Age group') +
  scale_fill_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12)



#####################################################
### Plot epi outputs
#####################################################
monthly_epi <- all_summ$monthly_epi_output_summ

# long datasets
prev_monthly <- monthly %>%
  dplyr::select(time, contains('prev'), month, ur, country, parameter_draw,
                site_name, target_type, parameter_draw) %>%
  pivot_longer(cols = contains('prev'),
               names_to = c(".value", "age_group_prev"),
               names_pattern = "(lm_prevalence)_(.*)") %>%
  distinct %>%
  mutate(date = as.Date(format(floor_date(date_decimal(time)), '%Y-%m-%d')),
         date = floor_date(date, unit ='month')) %>%
  filter(age_group_prev %in% c('5_9','9_18')) %>%
  mutate(age_group = case_when(
           age_group_prev == '5-9' ~ '5-8',
           age_group_prev == '9_18' ~ '9-17',
           TRUE ~ NA
         )) %>%
  filter(time >= 2025)

ggplot(prev_monthly) +
  geom_line(aes(x = date, y = lm_prevalence, group = age_group_prev, color = age_group_prev)) +
  geom_pointrange(data = pfpr %>% filter(age_group == '9-17' | age_group == '5-8'),
                  aes(x = date,
                      y = pf_positivity_rate/100,
                      ymin = pf_positivity_rate_lower/100,
                      ymax = pf_positivity_rate_upper/100,
                      group = age_group,
                      color = age_group),
                  position = position_dodge(width = 20),
                  size = 0.3) +
  facet_wrap(vars(country), scales = 'free')+
  scale_color_manual(values = ltc_cols_type) +
  scale_x_date(labels = scales::label_date_short()) +
  labs(x = 'Time',
       y = 'LM PfPR',
       color = NULL) +
  theme_classic(base_size = 12)
