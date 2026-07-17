
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
  select(country, site_name, ur,,
         time, infectivity_under5, infectivity_SAC, infectivity_16plus,
         mean_inf_under5, mean_inf_SAC, mean_inf_16plus,
         prop_mean_inf_under5, prop_mean_inf_SAC, prop_mean_inf_16plus,
         prop_sum_inf_under5, prop_sum_inf_SAC, prop_sum_inf_16plus
  ) %>%
  pivot_longer(
    cols = -c(time, site_name, country, ur),
    names_to = c(".value", "age_group"),
    names_pattern = "(infectivity|sum_inf|prop_sum_inf|mean_inf|prop_mean_inf)_(under5|SAC|16plus)"
  )

p1 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = infectivity, color = age_group)) +
  labs(y = 'Infectivity sum by age group',
       x = 'Year',
       color=  'Age group') +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

p2 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = mean_inf, color = age_group)) +
  labs(y = 'Mean infectivity per person',
       x = 'Year',
       color=  'Age group') +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

p3 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = prop_sum_inf, color = age_group)) +
  labs(y = 'Proportion of sum infectivity by age group',
       x = 'Year',
       color=  'Age group') +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

p4 <- ggplot(inf_long) +
  geom_line(aes(x = time, y = prop_mean_inf, color = age_group)) +
  labs(y = 'Proportion of mean infectivity by age group',
       x = 'Year',
       color=  'Age group') +
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))


# Summarize infectivity over last year of sim
infectivity_summ <- infectivity_annual_all %>%
  mutate(time = year) %>%
  filter(time == max(inf$time)) %>%
  select(country, site_name, ur,
         time, infectivity_under5, infectivity_SAC, infectivity_16plus,
         mean_inf_under5, mean_inf_SAC, mean_inf_16plus,
         prop_mean_inf_under5, prop_mean_inf_SAC, prop_mean_inf_16plus,
         prop_sum_inf_under5, prop_sum_inf_SAC, prop_sum_inf_16plus) %>%
  pivot_longer(
    cols = -c(time, site_name, country, ur),
    names_to = c(".value", "age_group"),
    names_pattern = "(infectivity|sum_inf|prop_sum_inf|mean_inf|prop_mean_inf)_(under5|SAC|16plus)"
  )

p5 <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_sum_inf), fill = 'darkred') +
  geom_text(aes(x = age_group, y = prop_sum_inf + 0.04, label = round(prop_sum_inf,2))) +
  labs(y = 'Proportion of sum infectivity per person by age group - not pop weighteds',
       x = 'Age group',
       subtitle = 'in last year of sim')+
  theme_classic(base_size = 12) +
  facet_wrap(vars(country))

p6 <- ggplot(infectivity_summ) +
  geom_col(aes(x = age_group, y = prop_mean_inf), fill = 'darkred') +
  geom_text(aes(x = age_group, y = prop_mean_inf + 0.04, label = round(prop_mean_inf,2))) +
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
print(p4)
print(p5)
print(p6)

# Close the PDF device to finalize the file
dev.off()


# Cases and prevalence
p7 <- ggplot(annual) +
  geom_line(aes(x = time, y = clinical, group = age_group, color = age_group)) +
  facet_wrap(vars(country))+
  labs(y = 'Clinical incidence per person per day, averaged per year',
       x = 'Time',
       color = 'Age group')+
  theme_classic(base_size = 12)

prev <- annual %>%
  select(time, contains('prev'), year, ur, country, site_name, parameter_draw) %>%
  pivot_longer(cols = contains('prev'),
               names_to = c(".value", "age_group_prev"),
               names_pattern = "(pcr_prevalence)_(0_5|2_10|5_16|16_100|0_100)") %>%
  distinct()

p8 <- ggplot(prev %>% filter(age_group_prev =='2_10')) +
  geom_line(aes(x = time, y = pcr_prevalence, group = age_group_prev, color = age_group_prev)) +
  facet_wrap(vars(country)) +
  labs(x = 'Time',
       y = 'PCR prevalence',
       color = 'Age group') +
  theme_classic(base_size = 12)

pdf(file = "outputs/epi_all_sites_annual2.pdf")

# Generate plots
print(p7)
print(p8)


# Close the PDF device to finalize the file
dev.off()
