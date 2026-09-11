# Data processing of Jen Hume's data on parasitaemia and entomology
# from August 2026

# Read in data
files <- list.files('data/Ento-reports-Aug2026', full.names = TRUE)
files <- files[grepl('Entomology', files)]
dfs <- lapply(files, function(x){
  d <- readr::read_csv(x)
  d$country <- str_extract(basename(x), "(?<=Entomology\\.Report\\.)[^.]+")
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
  mutate(month = str_sub(month_year, 1, 3),
         year = as.numeric(str_sub(month_year, 5, 8)),
         date = my(paste(month, year, sep = '-'))) %>%
  select(-month_year, -day_range) %>%
  # Group by month  (at the moment, first and second half of month are in separate lines)
  dplyr::summarise(
    dplyr::across(starts_with('n_'),
                  .fns = \(x) sum(x, na.rm = TRUE)),
    .by = c('month', 'year', 'date', 'country', 'age_group')
  ) %>%
  # Remake the rate variables and add in binomial CIs
  rowwise() %>%
  filter(n_dsf_performed > 0) %>%
  mutate(n_fed = n_mosq_in_cups - n_unfed,

         feeding_rate = round(n_fed / n_mosq_in_cups * 100, 1),
         feeding_rate_lower = binom.test(n_fed, n_mosq_in_cups)$conf.int[1] * 100,
         feeding_rate_upper = binom.test(n_fed, n_mosq_in_cups)$conf.int[2] * 100,

         survival_rate = round(n_mosq_dissected / n_fed * 100, 1),
         survival_rate_lower = binom.test(n_mosq_dissected, n_fed)$conf.int[1] * 100,
         survival_rate_upper = binom.test(n_mosq_dissected, n_fed)$conf.int[2] * 100,

         dsf_positivity_rate = round(n_dsf_positive / n_dsf_performed * 100, 1),
         dsf_positivity_rate_lower = binom.test(n_dsf_positive, n_dsf_performed)$conf.int[1] * 100,
         dsf_positivity_rate_upper = binom.test(n_dsf_positive, n_dsf_performed)$conf.int[2] * 100,

         mosq_positivity_rate = round(n_positive_mosq / n_mosq_dissected * 100, 1),
         mosq_positivity_rate_lower = binom.test(n_positive_mosq, n_mosq_dissected)$conf.int[1] * 100,
         mosq_positivity_rate_upper = binom.test(n_positive_mosq, n_mosq_dissected)$conf.int[2] * 100) %>%
  ungroup()


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
    facet_wrap(~country, scales = 'free') +
    labs(x = 'Date',
         y = paste(y_label, " (%)"),
         color = 'Age group') +
    theme_bw(base_size = 14)
}

# Usage - just provide the rate variable
plot_rate(dfs_df, "feeding_rate")
plot_rate(dfs_df, "survival_rate")
plot_rate(dfs_df, "dsf_positivity_rate")
plot_rate(dfs_df, "mosq_positivity_rate")

# Save data
saveRDS(dfs_df, 'data/Ento-reports-Aug2026/entomology_summarized.rds')

