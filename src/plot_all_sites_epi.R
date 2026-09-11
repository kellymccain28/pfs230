# Script to plot incidence and prevalence from all sites
library(tidyverse)

# Read in data from Jen
pfpr <- read_rds('data/Ento-reports-Aug2026/parasitaemia_summarized.rds')

# Plotting all sites together
path <- 'outputs/2026-09-07/' # first round of new calibration with Jen's data
all <- readRDS(paste0(path, 'all_processed_output.rds'))

outputs <- all %>%
  map('raw_output') %>%
  list_rbind()

daily <- all %>%
  map('daily_epi_output') %>%
  list_rbind() %>%
  mutate(age_upper = case_when(
    age_upper == 8.9972602739726 ~ 8,
    age_upper == 17.9972602739726 ~ 17,
    TRUE ~ age_upper),
  age_group = paste0(age_lower, '-', age_upper))

annual <- all %>%
  map('annual_epi_output') %>%
  list_rbind() %>%
  mutate(age_upper = case_when(
    age_upper == 8.9972602739726 ~ 8,
    age_upper == 17.9972602739726 ~ 17,
    TRUE ~ age_upper),
    age_group = paste0(age_lower, '-', age_upper))

monthly <- all %>%
  map('monthly_epi_output') %>%
  list_rbind() %>%
  mutate(age_upper = case_when(
    age_upper == 8.9972602739726 ~ 8,
    age_upper == 17.9972602739726 ~ 17,
    TRUE ~ age_upper),
    age_group = paste0(age_lower, '-', age_upper))

# long datasets
prev_annual <- annual %>%
  select(time, contains('prev'), year, ur, country,
         site_name, target_type, parameter_draw) %>%
  pivot_longer(cols = contains('prev'),
               names_to = c(".value", "age_group_prev"),
               names_pattern = "(lm_prevalence)_(.*)") %>%
  distinct()

prev_monthly <- monthly %>%
  dplyr::select(time, contains('prev'), month, ur, country,
         site_name, target_type, parameter_draw) %>%
  pivot_longer(cols = contains('prev'),
               names_to = c(".value", "age_group_prev"),
               names_pattern = "(lm_prevalence)_(.*)") %>%
  distinct %>%
  mutate(date = as.Date(format(floor_date(date_decimal(time)), '%Y-%m-%d')))

ltc_cols_age <- ltc::palettes$casa_natal
ltc_cols_type <- ltc::palettes$expevo

# Cases and prevalence
p7 <- ggplot(annual %>% filter(age_group %in% c('0-5','5-16','16-100')) %>%
               mutate(age_group = factor(age_group, levels = c('0-5','5-16','16-100')))) +
  geom_line(aes(x = time, y = clinical, group = age_group, color = age_group)) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))+
  scale_fill_manual(values = ltc_cols_age) +
  scale_color_manual(values = ltc_cols_age) +
  labs(y = 'Clinical incidence per person per day, averaged by year',
       x = 'Time',
       color = 'Age group',
       fill = 'Age group')+
  theme_classic(base_size = 12)


p8 <- ggplot(prev_annual %>% filter(age_group_prev =='2_10')) +
  geom_line(aes(x = time, y = lm_prevalence, group = target_type, color = target_type)) +
  facet_wrap(vars(country))+
  scale_color_manual(values = ltc_cols_type) +
  labs(x = 'Time',
       y = 'LM PfPR 2-10',
       color = NULL) +
  theme_classic(base_size = 12)

p9 <- ggplot(prev_monthly %>% filter(age_group_prev =='2_10')) +
  geom_line(aes(x = time, y = lm_prevalence, group = target_type, color = target_type)) +
  facet_wrap(vars(country))+
  scale_color_manual(values = ltc_cols_type) +
  labs(x = 'Time',
       y = 'LM PfPR 2-10',
       color = NULL) +
  theme_classic(base_size = 12)

# Make plots with new data
p10 <- prev_monthly %>%
  filter(age_group_prev =='5_9' & time > 2025) %>%
  ggplot() +
  geom_line(aes(x = date, y = lm_prevalence, group = target_type, color = target_type)) +
  geom_pointrange(data = pfpr %>% filter(age_group == '5-8'),
                  aes(x = date,
                      y = pf_positivity_rate/100,
                      ymin = pf_positivity_rate_lower/100,
                      ymax = pf_positivity_rate_upper/100,
                      group = age_group,
                      color = age_group),
                  position = position_dodge(width = 20),
                  size = 0.3) +
  facet_wrap(vars(country))+
  scale_color_manual(values = ltc_cols_type) +
  scale_x_date(labels = scales::label_date_short()) +
  labs(x = 'Time',
       y = 'LM PfPR 5-8',
       color = NULL) +
  theme_classic(base_size = 12)
ggsave(paste0(path, 'pfpr_5_8.png'), p10)

p11 <- prev_monthly %>%
  filter(age_group_prev =='9_18' & time > 2025) %>%
  ggplot() +
  geom_line(aes(x = date, y = lm_prevalence, group = target_type, color = target_type)) +
  geom_pointrange(data = pfpr %>% filter(age_group == '9-17'),
                  aes(x = date,
                      y = pf_positivity_rate/100,
                      ymin = pf_positivity_rate_lower/100,
                      ymax = pf_positivity_rate_upper/100,
                      group = age_group,
                      color = age_group),
                  position = position_dodge(width = 20),
                  size = 0.3) +
  facet_wrap(vars(country))+
  scale_color_manual(values = ltc_cols_type) +
  scale_x_date(labels = scales::label_date_short()) +
  labs(x = 'Time',
       y = 'LM PfPR 9-17',
       color = NULL) +
  theme_classic(base_size = 12)
ggsave(paste0(path, 'pfpr_9_17.png'), p11)

p12 <- prev_monthly %>%
  filter((age_group_prev =='9_18' | age_group_prev == '5_9') & time > 2025) %>%
  mutate(age_group_prev = ifelse(age_group_prev == '9_18', '9-17',
                                 ifelse(age_group_prev == '5_9', '5-8', NA))) %>%
  ggplot() +
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
  facet_wrap(vars(country))+
  scale_color_manual(values = ltc_cols_type[3:4]) +
  scale_x_date(labels = scales::label_date_short()) +
  labs(x = 'Time',
       y = 'LM PfPR',
       color = NULL) +
  theme_classic(base_size = 12)

pdf(file = paste0(path, "epi_all_sites_annual.pdf"), width = 11)

# Generate plots
print(p7)
print(p8)
print(p9)
print(p10)
print(p11)
print(p12)

# Close the PDF device to finalize the file
dev.off()

ggsave(paste0(path, 'lmprev2_10_allsites.png'), p8, width = 10)
ggsave(paste0(path, 'lmprev_allsites.png'), p12, width = 12, height = 5)
