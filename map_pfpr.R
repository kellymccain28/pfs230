# Combine MAP prevalence range outputs
# Manually chosen based on pixel map for each admin1 unit
# https://data.malariaatlas.org/trends?year=2024&metricGroup=Malaria&geographicLevel=pixel&metricSubcategory=Pf&metricType=rate&metricName=PR
library(janitor)
library(tidyverse)

files <- list.files('./site_files/map_pfpr_ranges/site-specific/', full.names = TRUE)
admin1 <- read.csv('./site_files/map_pfpr_ranges/map_pfpr_admin1.csv') %>%
  clean_names() %>%
  rename(site_name = name) %>%
  mutate(range = 'central',
         value = value / 100) %>%
  select(-iso3, -national_unit, -admin_level, -units)

pfprs <- lapply(files, function(fname){
  d <- read.csv(fname)
  site_name <- gsub("_pfpr_range\\.csv$", "", basename(fname[1]))
  if(site_name == 'greater-accra') site_name <- 'greater accra'
  d$site_name <- str_to_title(site_name)
  d
}) %>%
  bind_rows() %>%
  clean_names() %>%
  filter(metric == 'Infection Prevalence') %>%
  group_by(site_name) %>%
  mutate(
    lat_lowest  = latitude[year == 2024][which.min(value[year == 2024])],
    lat_highest = latitude[year == 2024][which.max(value[year == 2024])],
    range = case_when(
      latitude == lat_lowest  ~ 'lower',
      latitude == lat_highest ~ 'upper',
      TRUE ~ NA_character_
    ),
    value = value / 100
  ) %>%
  select(-lat_lowest, -lat_highest, -latitude, -longitude) %>%
  ungroup()

pfprs <- rbind(pfprs, admin1)

saveRDS(pfprs, './site_files/map_pfpr_ranges/map_pfpr_ranges.rds')

ggplot(pfprs %>%
         select(site_name, year, range, value) %>%
         pivot_wider(names_from = range, values_from = value)) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, color = site_name, fill = site_name),
              linetype = 2,
              alpha = 0.3) +
  geom_line(aes(x = year, y = central, color = site_name)) +
  geom_point(aes(x = year, y = central, color = site_name)) +
  facet_wrap(~site_name) +
  labs(title = 'MAP range of 2024 prevalence within admin 1 units',
       y = 'PfPR',
       fill = 'Site name',
       color = 'Site name',
       x = 'Year')

