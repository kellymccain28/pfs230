combine_summarize_dfs <- function(path){
  # With output, need to collate all of each type and save them
  # path <- 'outputs/2026-09-08/model_outputs/'
  files <- list.files(path, pattern = "^processed_out_.*\\.rds$", full.names = TRUE)

  # helper to pull one element from every file and rbind
  combine_element <- function(files, element) {
    map_dfr(files, function(f) {
      x <- readRDS(f)
      df <- x[[element]]
      df
    })
  }

  annual_all    <- combine_element(files, "annual_epi_output")
  monthly_all   <- combine_element(files, "monthly_epi_output")
  daily_all     <- combine_element(files, "daily_epi_output")
  infect_annual <- combine_element(files, "infectivity_annual")
  infect_month  <- combine_element(files, "infectivity_monthly")
  infect_all    <- combine_element(files, "infectivity")

  # one-row-per-draw parameter table, collected from the per-draw model_input
  params_all <- lapply(files, function(x){
    d <- readRDS(x)$model_input
    })
  names(params_all) <- sapply(params_all, function(d){
    paste(d$country, d$site_name, d$ur, d$target_type, d$parameter_draw, sep = '_')
  })

  # save the combined, tidy versions — this is what you'll actually load for analysis
  saveRDS(annual_all,    paste0(path, "combined_annual_epi_output.rds"))
  saveRDS(monthly_all,   paste0(path, "combined_monthly_epi_output.rds"))
  saveRDS(daily_all,     paste0(path, "combined_daily_epi_output.rds"))
  saveRDS(infect_annual, paste0(path, "combined_infectivity_annual.rds"))
  saveRDS(infect_month,  paste0(path, "combined_infectivity_monthly.rds"))
  saveRDS(infect_all,    paste0(path, "combined_infectivity.rds"))
  saveRDS(params_all,    paste0(path, "combined_params.rds"))



  # Then with those collated datasets, summarize over all parameter draws
  summarize_epi_df <- function(x){
    x %>%
      group_by(year, age_upper, age_lower, month, time,
               country, ur, site_name, burnin, target_type, week, day) %>%
      summarize(across(c(starts_with('lm_p'),
                         clinical, severe, mortality, yll, yld, dalys),
                       .fns = list(
                         median = ~median(.x, na.rm = TRUE),
                         lower  = ~quantile(.x, 0.025, na.rm = TRUE),
                         upper  = ~quantile(.x, 0.975, na.rm = TRUE)
                       ),
                       .names = "{.col}_{.fn}"),
                .groups = 'drop')
  }

  summarize_infec_df <- function(x){
    x %>%
      group_by(date, month, year,
               ur, country, site_name, burnin, target_type) %>%
      summarize(across(c(starts_with('infectivity'),
                         starts_with('n_age'), starts_with('mean_inf'),
                         starts_with('prop_')),
                       .fns = list(
                         median = ~median(.x, na.rm = TRUE),
                         lower  = ~quantile(.x, 0.025, na.rm = TRUE),
                         upper  = ~quantile(.x, 0.975, na.rm = TRUE)
                       ),
                       .names = "{.col}_{.fn}"),
                .groups = 'drop')
  }

  annual_all_summ <- annual_all %>%
    summarize_epi_df()

  monthly_all_summ <- monthly_all %>%
    summarize_epi_df()

  daily_all_summ <- daily_all %>%
    summarize_epi_df()

  infect_annual_summ <- infect_annual %>%
    summarize_infec_df()

  infect_monthly_summ <- infect_month %>%
    summarize_infec_df()

  infect_daily_summ <- infect_all %>%
    summarize_infec_df()

  # save the combined, tidy versions — this is what you'll actually load for analysis
  saveRDS(annual_all_summ,    paste0(path, "summarized_annual_epi_output.rds"))
  saveRDS(monthly_all_summ,   paste0(path, "summarized_monthly_epi_output.rds"))
  saveRDS(daily_all_summ,     paste0(path, "summarized_daily_epi_output.rds"))
  saveRDS(infect_annual_summ, paste0(path, "summarized_infectivity_annual.rds"))
  saveRDS(infect_monthly_summ,  paste0(path, "summarized_infectivity_monthly.rds"))
  saveRDS(infect_daily_summ,    paste0(path, "summarized_infectivity.rds"))

  processed_output_all <- list('annual_epi_output_summ' = annual_all_summ,
                               'monthly_epi_output_summ' = monthly_all_summ,
                               'daily_epi_output_summ' = daily_all_summ,
                               'infectivity_annual_summ' = infect_annual_summ,
                               'infectivity_monthly_summ' = infect_monthly_summ,
                               'infectivity_daily_summ' = infect_daily_summ)

  saveRDS(processed_output_all, paste0(path, "processed_output_summ.rds"))
}
