
# Plotting all sites together
all <- readRDS('M:/Kelly/postdoc_JoeC/pfs230/outputs/all_processed_output.rds')
outputs <- all %>%
  map('raw_output') %>%
  list_rbind()

infectivity_annual_all <- all %>%
  map('infectivity_annual') %>%
  list_rbind()

infectivity_daily_all <- all %>%
  map('infectivity') %>%
  list_rbind()

daily <- all %>%
  map('daily_epi_output') %>%
  list_rbind() %>%
  mutate(age_group = paste0(age_lower, '-', age_upper))

annual <- all %>%
  map('annual_epi_output') %>%
  list_rbind() %>%
  mutate(age_group = paste0(age_lower, '-', age_upper))

ltc_cols_age <- ltc::palettes$casa_natal
ltc_cols_type <- ltc::palettes$expevo

# Get long df for annual data
inf_long <- infectivity_annual_all %>%
  mutate(time = year) %>%
  select(country, site_name, ur, target_type,
         time, infectivity_under5, infectivity_SAC, infectivity_16plus,
         mean_inf_under5, mean_inf_SAC, mean_inf_16plus,
         prop_mean_inf_under5, prop_mean_inf_SAC, prop_mean_inf_16plus,
         prop_sum_inf_under5, prop_sum_inf_SAC, prop_sum_inf_16plus
  ) %>%
  pivot_longer(
    cols = -c(time, site_name, country, ur, target_type),
    names_to = c(".value", "age_group"),
    names_pattern = "(infectivity|sum_inf|prop_sum_inf|mean_inf|prop_mean_inf)_(under5|SAC|16plus)"
  ) %>%
  mutate(age_group = factor(age_group, levels = c('under5','SAC','16plus')),
         target_type = factor(target_type, levels = c('lower','central','upper')))

# Get long df for daily data
inf_long_daily <- infectivity_daily_all %>%
  select(country, site_name, ur, target_type,
         timestep, year, infectivity_under5, infectivity_SAC, infectivity_16plus,
         mean_inf_under5, mean_inf_SAC, mean_inf_16plus,
         prop_mean_inf_under5, prop_mean_inf_SAC, prop_mean_inf_16plus,
         prop_sum_inf_under5, prop_sum_inf_SAC, prop_sum_inf_16plus
  ) %>%
  pivot_longer(
    cols = -c(timestep, year, site_name, country, ur, target_type),
    names_to = c(".value", "age_group"),
    names_pattern = "(infectivity|sum_inf|prop_sum_inf|mean_inf|prop_mean_inf)_(under5|SAC|16plus)"
  ) %>%
  mutate(age_group = factor(age_group, levels = c('under5','SAC','16plus')),
         target_type = factor(target_type, levels = c('lower','central','upper')))

# this is the total age-specific daily infectivity averaged per year
p1 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = infectivity, color = target_type, linetype = target_type)) +
  labs(y = 'Total infectivity by age group',
       x = 'Year',
       color=  'Target type',
       linetype = 'Target type') +
  scale_color_manual(values = ltc_cols_type) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(age_group))

# mean infectivity per person per year
p2 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = mean_inf, color = age_group, linetype = target_type)) +
  labs(y = 'Mean infectivity per person',
       x = 'Year',
       color = 'Age group',
       linetype = NULL) +
  scale_color_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(age_group))

# Infectivity daily for last 5 years of sim
inf_recent <- inf_long_daily %>%
  filter(year >= (max(inf_long_daily$year)-5))

p1b <- ggplot(inf_recent) +
  geom_line(aes(x = timestep/365 + 2000, y = infectivity, color = age_group)) +
  labs(y = 'Total daily infectivity by age group',
       x = 'Year',
       color=  'Age group') +
  scale_color_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

p2b <- ggplot(inf_recent) +
  geom_line(aes(x = timestep/365 + 2000, y = mean_inf, color = age_group)) +
  labs(y = 'Mean infectivity per person',
       x = 'Year',
       color=  'Age group') +
  scale_color_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

# Annual *
p3 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = prop_sum_inf, color = age_group, linetype = target_type)) +
  labs(y = 'Proportion of total infectivity',
       x = 'Year',
       color=  'Age group') +
  scale_color_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

# *
p3b <- ggplot(inf_long) +
  geom_col(aes(x = time, y = prop_sum_inf, fill = age_group, color = age_group),
           position = position_fill(),
           width = 1) +
  scale_x_continuous(breaks = seq(2000,2030,5)) +
  labs(y = 'Proportion of total infectivity',
       x = 'Year',
       color =  'Age group',
       fill =  'Age group',
       caption = paste0()) +
  scale_color_manual(values = ltc_cols_age) +
  scale_fill_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

p4 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = prop_mean_inf, color = age_group)) +
  labs(y = 'Proportion of mean infectivity by age group',
       x = 'Year',
       color=  'Age group') +
  scale_color_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))


# Summarize infectivity over last year of sim
infectivity_lastyear_tbl <- inf_long %>%
  filter(time == 2026) %>%
  select(country, site_name, target_type, age_group, prop_sum_inf) %>%
  mutate(prop_sum_inf = round(prop_sum_inf, 2)) %>%
  pivot_wider(id_cols = c(country, site_name, target_type),
              names_from = age_group,
              values_from = prop_sum_inf) %>%
  mutate(target_type = factor(target_type, levels = c('lower','central','upper')))
write.csv(infectivity_lastyear_tbl, "outputs/prop_sum_infectivity.csv", row.names = FALSE)


infectivity_summ <- infectivity_annual_all %>%
  mutate(time = year) %>%
  filter(time == 2026) %>%
  select(country, site_name, ur,target_type,
         time, infectivity_under5, infectivity_SAC, infectivity_16plus,
         mean_inf_under5, mean_inf_SAC, mean_inf_16plus,
         prop_mean_inf_under5, prop_mean_inf_SAC, prop_mean_inf_16plus,
         prop_sum_inf_under5, prop_sum_inf_SAC, prop_sum_inf_16plus) %>%
  pivot_longer(
    cols = -c(time, site_name, country, ur, target_type),
    names_to = c(".value", "age_group"),
    names_pattern = "(infectivity|sum_inf|prop_sum_inf|mean_inf|prop_mean_inf)_(under5|SAC|16plus)"
  ) %>%
  mutate(age_group = factor(age_group, levels = c('under5','SAC','16plus')),
         target_type = factor(target_type, levels = c('lower','central','upper')))

# *
p5 <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_sum_inf, fill = age_group)) +
  geom_text(aes(x = age_group, y = prop_sum_inf + 0.04, label = round(prop_sum_inf,2)),
            size = 3) +
  labs(y = 'Proportion of sum infectivity',
       x = 'Age group',
       color = 'Age group')+
  scale_fill_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12)  +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

# *
p5b <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_sum_inf, fill = target_type),
           position = 'dodge') +
  geom_text(aes(x = age_group, y = prop_sum_inf + 0.03,
                group = target_type, label = round(prop_sum_inf,2)),
            position = position_dodge(width = .9), size = 3) +
  labs(y = 'Proportion of total infectivity',
       x = 'Age group',
       fill = NULL)+
  scale_fill_manual(values = ltc_cols_type) +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

p6a <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_mean_inf, fill = age_group)) +
  geom_text(aes(x = age_group, y = prop_mean_inf + 0.04, label = round(prop_mean_inf,2)),
            size = 3) +
  labs(y = 'Relative per-person infectivity',
       x = 'Age group',
       color = 'Age group')+
  scale_fill_manual(values = ltc_cols_age) +
  theme_classic(base_size = 12)  +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

p6 <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_mean_inf, fill = target_type),
           position = 'dodge') +
  geom_text(aes(x = age_group, y = prop_mean_inf + 0.034,
                group = target_type, label = round(prop_mean_inf,2)),
            position = position_dodge(width = .9), size = 3) +
  labs(y = 'Relative per-person infectivity',
       x = 'Age group',
       fill = NULL)+
  scale_fill_manual(values = ltc_cols_type) +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

p6c <- ggplot(infectivity_summ) +
  geom_col(aes(x = target_type, y = prop_mean_inf, fill = age_group),
           position = 'fill') +
  geom_text(aes(x = target_type, y = prop_mean_inf,
                group = age_group, label = round(prop_mean_inf,2)),
            position = position_fill(), size = 3) +
  labs(y = 'Relative per-person infectivity',
       x = 'Age group',
       fill = NULL)+
  scale_fill_manual(values = ltc_cols_type) +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

# Absolute values of mean and summed infectivity
p7 <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = infectivity, fill = target_type),
           position = 'dodge') +
  geom_text(aes(x = age_group, y = infectivity + 5,
                group = target_type, label = round(infectivity,1)),
            position = position_dodge(width = .9), size = 3) +
  labs(y = 'Total infectivity per age group',
       x = 'Age group',
       fill = NULL)+
  scale_fill_manual(values = ltc_cols_type) +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

p8 <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = mean_inf, fill = target_type),
           position = 'dodge') +
  geom_text(aes(x = age_group, y = mean_inf + 0.002,
                group = target_type, label = round(mean_inf,3)),
            position = position_dodge(width = .9), size = 3) +
  labs(y = 'Mean per-person infectivity',
       x = 'Age group',
       fill = NULL)+
  scale_fill_manual(values = ltc_cols_type) +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

# Save individual plots
ggsave('outputs/proportion_total_infectivity_typefacet.png', p5, width = 12)
ggsave('outputs/proportion_total_infectivity.png', p5b, width = 12)
ggsave('outputs/relative_per_person_infectivity_typefacet.png', p6a, width = 12)
ggsave('outputs/relative_per_person_infectivity.png', p6, width = 12)
ggsave('outputs/relative_per_person_infectivity_fill.png', p6c, width = 12)
ggsave('outputs/total_infectivity.png', p7, width = 12)
ggsave('outputs/mean_per_person_infectivity.png', p8, width = 12)


# Save all infectivity plots
pdf(file = "outputs/infectivity_all_sites_annual.pdf", width = 12)

# Generate plots
print(p1)
print(p2)
print(p1b)
print(p2b)
print(p3)
print(p3b)
print(p4)
print(p5)
print(p5b)
print(p6a)
print(p6)
print(p6c)
print(p7)
print(p8)

# Close the PDF device to finalize the file
dev.off()

