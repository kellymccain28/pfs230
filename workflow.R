library(hipercow)

#Prep hipercow:
hipercow_environment_create(sources = c("helper_functions.R",
                                        "create_site_file.R"))
hipercow_provision(method = 'script')


# Calibration to high and low values of parasite prevalence for each site
pr_match(site_name = 'Atlantique',
         pfpr_target_type = 'lower')
#test
PRmatch_draws_Atlantique_lower <- readRDS("M:/Kelly/postdoc_JoeC/pfs230/PrEIR/PRmatch_draws_Atlantique_lower.rds")
parameters <- set_equilibrium(params, init_EIR = PRmatch_draws_Atlantique_lower$starting_EIR)
raw <- run_simulation(parameters$timesteps + 100, parameters = parameters)
raw$pfpr <- raw$n_detect_lm_730_3650  / raw$n_age_730_3650
raw$year <- ceiling(raw$timestep / 365) + 2000 - 16 # assuming 15 years burnin
years <- unique(raw[raw$year %in% seq(2010,2024),]$year)
ggplot() +
  geom_point(aes(x = years, y = target/100), col = "dodgerblue", size = 4) +
  geom_line(data = raw, aes(x = timestep/365 + 1985, y = pfpr), col = "deeppink", linewidth = 1) +
  ylim(0, 1) +
  ylab(expression(italic(Pf)*Pr[2-10])) +
  xlab("Time") +
  theme_bw()
# need to make a function to do this for all of them at once


# Run simulations for all sites
t1 <- task_create_expr(expr = create_site_file(),
                 resources = hipercow_resources(cores = 6))
task_log_show(t1)
