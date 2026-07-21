# Script to calibrate and plot calibration for each site
source('helper_functions.R')

# Calibration to high and low values of parasite prevalence for each site
target_types <- c('upper', 'lower', 'central')
admin1s <- c('Centre-Sud','Greater Accra','Kisumu','Koulikoro','Pwani','Atlantique')

for(a in admin1s){
  for(t in target_types){
    pr_match(site_name = a,
             pfpr_target_type = t)
  }
}

#test the calibration
# files <- list.files('M:/Kelly/postdoc_JoeC/pfs230/PrEIR/', full.names = TRUE)
# preir <- bind_rows(lapply(files, readRDS))
# saveRDS(preir, 'M:/Kelly/postdoc_JoeC/pfs230/PrEIR/PRmatch_draws.rds')
preir <- readRDS('M:/Kelly/postdoc_JoeC/pfs230/PrEIR/PRmatch_draws.rds')
combos <- split(crossing(target_types, admin1s), ~ admin1s + target_types)
all_model_input <- readRDS("M:/Kelly/postdoc_JoeC/pfs230/site_files/all_model_input.rds")
map_pfpr_ranges <- readRDS("site_files/map_pfpr_ranges/map_pfpr_ranges.rds")
#site file outputs:
all_processed_output <- readRDS("M:/Kelly/postdoc_JoeC/pfs230/outputs/all_processed_output.rds")
daily <- all_processed_output %>%
  map('raw_output') %>%
  map(~ postie::get_prevalence(
    .x,
    diagnostic = 'lm',
    baseline_year = 2000
  ) %>%
    mutate(site_name = unique(.x$site_name))) %>%
  list_rbind()

# Compare site file outputs, cali calibrated to lower, upper, and central, and target MAP pfpr values
lapply(combos, function(cc){
  # Get site-specific calibrated EIRs to central, higher, or lower prevalence MAP estimates
  preir <- preir %>% filter(site_name == cc$admin1s & pfpr_target_type == cc$target_types)
  # Get MAP data and filter to only 2010-2024 (this is what I calibrated to to be consistent with {cali})
  map_data <- map_pfpr_ranges[map_pfpr_ranges$site_name == cc$admin1s & map_pfpr_ranges$range == cc$target_types,]
  target <- map_data[map_data$year %in% 2010:2024,]$value#c(seq(2000,2024,5),2024),]$value

  # Get parameters for specific site, using calibrated initial EIR
  params <- all_model_input[grep(cc$admin1s, names(all_model_input), value = TRUE, ignore.case = TRUE)][[1]]$param_list
  parameters <- set_equilibrium(params, init_EIR = preir$starting_EIR)

  message('running simulation for ', preir)
  raw <- run_simulation(parameters$timesteps, parameters = parameters)
  raw <- postie::drop_burnin(raw, burnin = 15 * 365)
  prev <- raw %>%
    postie::get_prevalence(
      diagnostic = 'lm',
      baseline_year = params$start_year + 15 # 2000
    )
  years <- unique(prev[prev$year %in% seq(2010,2024),]$year)

  p <- ggplot() +
    geom_point(aes(x = years, y = target, col = "MAP prevalence\ntarget"), size = 4) +
    geom_line(data = prev, aes(x = time, y = lm_prevalence_2_10, col = 'cali')) +
    ylim(0, 1) +
    labs(y = expression(italic(Pf)*Pr[2-10]),
         x = "Time",
         color = NULL) +
    theme_bw() +
    theme(legend.position = c(0.8,0.8))

  if(cc$target_types == 'central'){
    daily1 <- daily[daily$site_name == cc$admin1s,]
    p <- p +
      geom_line(data = daily1, aes(x = time, y = lm_prevalence_2_10, color = 'Site file'))
  }
  p
  ggsave(paste0('M:/Kelly/postdoc_JoeC/pfs230/outputs/preir_validation/plot_', cc$admin1s, '_', cc$target_types, '.png'), p)
}
)

# plot the raw model output using get_equilibrium() vs the central target (MAP admin1 level prevalence)
combos_central <- combos %>%
  bind_rows() %>% filter(target_types == 'central') %>%
  split(~ admin1s + target_types)

lapply(combos_central, function(cc){
  daily1 <- daily[daily$site_name == cc$admin1s,]
  map_data <- map_pfpr_ranges[map_pfpr_ranges$site_name == cc$admin1s & map_pfpr_ranges$range == cc$target_types,]
  target <- map_data[map_data$year %in% 2010:2024,]$value#c(seq(2000,2024,5),2024),]$value
  years = seq(2010,2024)

  p <- ggplot() +
    geom_point(aes(x = years, y = target), col = "dodgerblue", size = 4) +
    geom_line(data = daily1, aes(x = time, y = pcr_prevalence_2_10), col = "deeppink", linewidth = 1) +
    ylim(0, 1) +
    ylab(expression(italic(Pf)*Pr[2-10])) +
    xlab("Time") +
    theme_bw()

  ggsave(paste0('M:/Kelly/postdoc_JoeC/pfs230/outputs/preir_validation/plot_', cc$admin1s, '_', cc$target_types, '_sitefileEIR.png'), p)
})
