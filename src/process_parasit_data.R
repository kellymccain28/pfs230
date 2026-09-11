# Data processing of Jen Hume's data on parasitaemia and entomology
# from August 2026

# Read in data
files <- list.files('data/Ento-reports-Aug2026', full.names = TRUE)
files <- files[grepl('Parasit', files)]
dfs <- lapply(files, function(x){
  d <- readr::read_csv(x)
  d$country <- str_extract(basename(x), "(?<=Parasitology\\.Report\\.)[^.]+")
  return(d)
}
)

dfs_df <- bind_rows(dfs) %>%
  janitor::clean_names() %>%
  mutate(age_group = case_when(
    age_group == 'Total' ~ '5-17',
    TRUE ~ age_group)) %>%
  # align age groups with malsim output
  mutate(age_group = stringr::str_replace(age_group, 'yo','')) %>%
  # divide month and year into 2 variables
  mutate(month_lab = str_sub(month_year, 1, 3),
         year = as.numeric(str_sub(month_year, 5, 8)),
         date = my(paste(month_lab, year, sep = '-')),
         month = month(date)) %>%
  select(-month_year, -day_range) %>%
  # Group by month  (at the moment, first and second half of month are in separate lines)
  dplyr::summarise(
    dplyr::across(starts_with('n_'),
                  .fns = \(x) sum(x, na.rm = TRUE)),
    .by = c('month', 'year', 'date', 'country', 'age_group')
  ) %>%
  # Remake the rate variables and add in binomial CIs
  rowwise() %>%
  mutate(pf_positivity_rate = ifelse(n_volunteers > 0, round(n_pf_positive / n_volunteers * 100, 1), 0),
         pf_positivity_rate_lower = ifelse(n_volunteers > 0, binom.test(n_pf_positive, n_volunteers)$conf.int[1] * 100, 0),
         pf_positivity_rate_upper = ifelse(n_volunteers > 0, binom.test(n_pf_positive, n_volunteers)$conf.int[2] * 100, 0),

         gam_positivity_rate = ifelse(n_volunteers > 0, round(n_gam_positive / n_volunteers * 100, 1), 0),
         gam_positivity_rate_lower = ifelse(n_volunteers > 0, binom.test(n_gam_positive, n_volunteers)$conf.int[1] * 100, 0),
         gam_positivity_rate_upper = ifelse(n_volunteers > 0, binom.test(n_gam_positive, n_volunteers)$conf.int[2] * 100, 0)) %>%
  ungroup() %>%
  # Add site names (BF and Mali are not in this df, but included for completeness)
  mutate(site_name = case_when(
    country == 'Benin' ~ 'Atlantique',
    country == 'Burkina Faso' ~ 'Centre-Sud',
    country == 'Ghana' ~ 'Greater Accra',
    country == 'Kenya' ~ 'Kisumu',
    country == 'Mali' ~ 'Koulikoro',
    country == 'Tanzania' ~ 'Pwani',
    TRUE ~ NA
  ))


# Plot to check CIs
plot_rate <- function(data, rate_var) {

  lower_var <- paste0(rate_var, "_lower")
  upper_var <- paste0(rate_var, "_upper")

  y_label <- stringr::str_to_sentence(gsub("_", " ", rate_var))

  # Create the plot
  ggplot(data, aes(x = date,
                   y = .data[[rate_var]],
                   ymin = .data[[lower_var]],
                   ymax = .data[[upper_var]],
                   group = age_group,
                   color = age_group)) +
    geom_pointrange(position = position_dodge(width = 20),
                    size = 0.3) +
    facet_wrap(~country) +
    labs(x = 'Date',
         y = paste(y_label, " (%)"),
         color = 'Age group') +
    theme_bw(base_size = 14)
}

# Usage - just provide the rate variable
plot_rate(dfs_df, "pf_positivity_rate")
plot_rate(dfs_df, "gam_positivity_rate")

# Save data
saveRDS(dfs_df, 'data/Ento-reports-Aug2026/parasitaemia_summarized.rds')


