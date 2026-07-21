
# Plotting all sites together
all <- readRDS('M:/Kelly/postdoc_JoeC/pfs230/outputs/all_processed_output.rds')
outputs <- all %>%
  map('raw_output') %>%
  list_rbind()

infectivity_annual_all <- all %>%
  map('infectivity_annual') %>%
  list_rbind()

daily <- all %>%
  map('daily_epi_output') %>%
  list_rbind() %>%
  mutate(age_group = paste0(age_lower, '-', age_upper))

annual <- all %>%
  map('annual_epi_output') %>%
  list_rbind() %>%
  mutate(age_group = paste0(age_lower, '-', age_upper))



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
  )

p1 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = infectivity, color = age_group, linetype = target_type)) +
  labs(y = 'Infectivity sum by age group',
       x = 'Year',
       color=  'Age group') +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(age_group))

p2 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = mean_inf, color = age_group, linetype = target_type)) +
  labs(y = 'Mean infectivity per person',
       x = 'Year',
       color=  'Age group',
       linetype = NULL) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(age_group))

# *
p3 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = prop_sum_inf, color = age_group, linetype = target_type)) +
  labs(y = 'Proportion of sum infectivity by age group',
       x = 'Year',
       color=  'Age group') +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

# *
p3b <- ggplot(inf_long) +
  geom_col(aes(x = time, y = prop_sum_inf, fill = age_group, color = age_group),
           position = position_fill(),
           width = 1) +
  scale_x_continuous(breaks = seq(2000,2030,5)) +
  labs(y = 'Proportion of sum infectivity by age group',
       x = 'Year',
       color =  'Age group',
       fill =  'Age group',
       caption = paste0()) +
  theme_classic(base_size = 12) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

p4 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = prop_mean_inf, color = age_group)) +
  labs(y = 'Proportion of mean infectivity by age group',
       x = 'Year',
       color=  'Age group') +
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
              values_from = prop_sum_inf)
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
  )

# *
p5 <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_sum_inf), fill = 'darkred') +
  geom_text(aes(x = age_group, y = prop_sum_inf + 0.04, label = round(prop_sum_inf,2))) +
  labs(y = 'Proportion of sum infectivity per person by age group - not pop weighted',
       x = 'Age group',
       subtitle = 'in last year of sim')+
  theme_classic(base_size = 12)  +
  facet_grid(rows = vars(country),
             cols = vars(target_type))

# *
p5b <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_sum_inf, fill = target_type),
           position = 'dodge') +
  geom_text(aes(x = age_group, y = prop_sum_inf + 0.03,
                group = target_type, label = round(prop_sum_inf,2)),
            position = position_dodge(width = .9)) +
  labs(y = 'Proportion of sum infectivity per person by age group - not pop weighted',
       x = 'Age group',
       subtitle = 'in last year of sim',
       fill = NULL)+
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

p6 <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_mean_inf, fill = target_type),
           position = 'dodge') +
  geom_text(aes(x = age_group, y = prop_mean_inf + 0.034,
                group = target_type, label = round(prop_mean_inf,2)),
            position = position_dodge(width = .9)) +
  labs(y = 'Proportion of mean infectivity per person by age group - not pop weighteds',
       x = 'Age group',
       subtitle = 'in last year of sim')+
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

pdf(file = "outputs/infectivity_all_sites_annual.pdf")

# Generate plots
print(p1)
print(p2)
print(p3)
print(p3b)
print(p4)
print(p5)
print(p5b)
print(p6)

# Close the PDF device to finalize the file
dev.off()


# Cases and prevalence
p7 <- ggplot(annual) +
  geom_line(aes(x = time, y = clinical, group = age_group, color = age_group)) +
  facet_grid(rows = vars(country),
             cols = vars(target_type))+
  labs(y = 'Clinical incidence per person per day, averaged per year',
       x = 'Time',
       color = 'Age group')+
  theme_classic(base_size = 12)

prev <- annual %>%
  select(time, contains('prev'), year, ur, country,
         site_name, target_type, parameter_draw) %>%
  pivot_longer(cols = contains('prev'),
               names_to = c(".value", "age_group_prev"),
               names_pattern = "(lm_prevalence)_(0_5|2_10|5_16|16_100|0_100)") %>%
  distinct()

p8 <- ggplot(prev %>% filter(age_group_prev =='2_10')) +
  geom_line(aes(x = time, y = lm_prevalence, group = target_type, color = target_type)) +
  facet_wrap(vars(country))+
  labs(x = 'Time',
       y = 'LM prevalence',
       color = 'Age group') +
  theme_classic(base_size = 12)

pdf(file = "outputs/epi_all_sites_annual2.pdf")

# Generate plots
print(p7)
print(p8)


# Close the PDF device to finalize the file
dev.off()
