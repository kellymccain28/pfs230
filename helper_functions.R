# Functions to run the TBV model (based on JDC's function first_try_gather_function.R)

# Function to fetch site files and save them
# input: vector of country iso3 country codes and path to save?
#' @param country_code iso3c country code, scalar character
#' @param admin_1_name name of the admin 1 unit, scalar character
#' @param ur urban or rural, character
#' @param hum_pop human poplulation for model runs
#' @param quick_run logical; if T, will parameterise for smaller pop for shorter period
fetch_all_sites <- function(country_code,
                            admin_1_name,
                            ur
){

  sites <- site::fetch_site(iso3c = country_code,
                               admin_level = 1,
                               urban_rural = TRUE)

  filter_site_by_region <- function(sites, admin_1_name, ur = NULL) {
    if (is.data.frame(sites)) {
      if ("name_1" %in% colnames(sites)) {
        nm <- sites$name_1
        Encoding(nm) <- "UTF-8"
        nm <- iconv(nm, from = "UTF-8", to = "ASCII//TRANSLIT")
        keep <- nm == admin_1_name
        if (!is.null(ur) && "urban_rural" %in% colnames(sites)) {
          keep <- keep & sites$urban_rural == ur
        }
        sites <- sites[keep, , drop = FALSE]
      }
      return(sites)
    } else if (is.list(sites)) {
      return(lapply(sites, filter_site_by_region, admin_1_name = admin_1_name, ur = ur))
    } else {
      return(sites)
    }
  }
    # Filter to only specific site
  site <- filter_site_by_region(sites,
                                admin_1_name = admin_1_name,
                                ur = ur)


  #### Save the site file
  if (!dir.exists('site_files')) dir.create('site_files')
  fname <- paste0('site_files/site_', country_code, '_', admin_1_name, '_', ur, '.rds')
  saveRDS(site, fname)
  message('saved ', fname)

  return(site)
}


# Diagnostics for site files - to do
# function to look at the interventions over time
diagnostics_site_files <- function(site){

  # Plot interventions over time
}

# Function to determine the age groups and time horizons
pull_age_groups_time_horizon<- function(quick_run = T){

  year<- 365
  burnin<- 15

  if(quick_run == TRUE){

    term_yr<- 2026
    pop_val<- 5000

    min_ages = c(0, 2, 5, 16, 0)*year
    max_ages = c(5, 10, 16, 100, 100) * year

  } else{

    pop_val<- 50000
    term_yr<- 226

    min_ages = c(0, 2, 5, 16, 0)*year
    max_ages = c(5, 10, 16, 100, 100) * year

  }

  return(list('term_yr' = term_yr,
              'pop_val' = pop_val,
              'min_ages'= min_ages,
              'max_ages' = max_ages,
              'burnin' = burnin))
}

# Function to generate baseline parameters based on site files
gather_params <- function(site, # this will be the output of the fetch_all_sites() that had been saved
                          hum_pop = 20000,
                          quick_run = TRUE,
                          parameter_draw = 0
){


  run_params<- pull_age_groups_time_horizon(quick_run = quick_run)

  # Set up timing and magnitude of ITN distribution
  site$interventions$itn$implementation$itn_input_dist <- site::site_usage_to_model_distribution(
    usage = site$interventions$itn$use$itn_use,
    usage_year = site$interventions$itn$use$year,
    usage_day_of_year = site$interventions$itn$use$usage_day_of_year,
    distribution_year = site$interventions$itn$implementation$year,
    distribution_day_of_year = site$interventions$itn$implementation$distribution_day_of_year,
    distribution_lower = site$interventions$itn$implementation$distribution_lower,
    distribution_upper = site$interventions$itn$implementation$distribution_upper,
    net_loss_function = netz::net_loss_map,
    half_life = site$interventions$itn$retention_half_life)

  # pull parameters for this site ------------------------------------------------
  params <- site::site_parameters(
    interventions = site$interventions,
    demography = site$demography,
    vectors = site$vectors,
    seasonality = site$seasonality,
    eir = site$eir$eir[1],

    age_group = run_params$min_ages, # ?
    start_year = min(site$demography$year) - run_params$burnin,
    end_year = max(site$prevalence$year),
    overrides = list(human_population = run_params$pop_val)
  )


  # set age groups
  params$clinical_incidence_rendering_min_ages = run_params$min_ages
  params$clinical_incidence_rendering_max_ages = run_params$max_ages
  params$severe_incidence_rendering_min_ages = run_params$min_ages
  params$severe_incidence_rendering_max_ages = run_params$max_ages
  params$age_group_rendering_min_ages = run_params$min_ages
  params$age_group_rendering_max_ages = run_params$max_ages
  params$prevalence_rendering_min_ages = run_params$min_ages
  params$prevalence_rendering_max_ages = run_params$max_ages

  # if this is a stochastic run, set parameter draw ------------------------------
  if (parameter_draw > 0){
    params<- params |>
      malariasimulation::set_parameter_draw(parameter_draw)
  }

  #params$pev<- TRUE # ?
  inputs <- list( # Makes more sense inside a function
    'param_list' = params,
    'site_name' = site$sites$name_1,
    'ur' = site$sites$urban_rural,
    'country' = site$sites$country,
    #'scenario' = scenario,
    'parameter_draw' = parameter_draw,
    'pop_val' = run_params$pop_val,
    'burnin' =  run_params$burnin
  )
  return(inputs)
}

# Function to add in the vaccines to the intervention parameters

# Target function for calibration -- will do this later
# calibration to average annual pfpr value for last 5 years of simulation
# annual_pfpr_summary <- function(x){
#
#   x$year <- ceiling(x$timestep / 365)
#   x <- x[x$year == max(x$year) | x$year == max(x$year)-1,]
#   pfpr <- x$n_detect_730_3650 / x$n_730_3650
#   year <- x$year
#   tapply(pfpr, year, mean)
# }

# Calibration of the model to site-specific information -- to do later



# Function to run model for the individual
run_model<- function(model_input,
                     verbose = F){

  params <- model_input$param_list
  params$progress_bar <- TRUE

  if(verbose==T){
    print(params)
  }

  # Set equilibrium
  message('set equilibrium')
  params_equil <- malariasimulation::set_equilibrium(
    params,
    init_EIR = params$init_EIR
  )

  message('running the model')
  model <- retry::retry(
    malariasimulation::run_simulation(timesteps = params_equil$timesteps,
                                      parameters = params_equil),
    max_tries = 5,
    when = 'error reading from connection|embedded nul|unknown type',
    interval = 3
  )

  # add identifying information to output
  model <- model |>
    dplyr::mutate(country = model_input$country,
                  site_name = model_input$site_name,
                  urban_rural = model_input$ur,
                  population = model_input$pop_val,
                  parameter_draw = model_input$parameter_draw,
                  burnin = model_input$burnin)
  #iso = model_input$iso3c,
  #description = model_input$description,
  #scenario = model_input$scenario,
  #gfa = model_input$gfa,

  # save model runs somewhere # ?
  message('saving the model')
  return(model)
}

# Processing the model output from run_model()
# outputs a processed data frame
process_output <- function(model){

  # Drop burnin
  raw_output <- postie::drop_burnin(model, burnin = run_parameters$burnin)

  message('calculating rates')

  rates <- postie::get_rates(
    raw_output,
    baseline_year = 1
  ) %>%# add identifying information to output
    mutate(scen_name = scen_name_val,
           pfpr = pfpr_val,
           seas_name = seas_name_val,
           parameter_draw = run_parameters$parameter_draw,
           population = run_parameters$population,
           burnin = run_parameters$burnin)

  rates_annual <- rates %>%
    dplyr::summarise(
      clinical = stats::weighted.mean(clinical, person_days),
      severe = stats::weighted.mean(severe, person_days),
      mortality = stats::weighted.mean(mortality, person_days),
      yll = stats::weighted.mean(yll, person_days),
      yld = stats::weighted.mean(yld, person_days),
      dalys = stats::weighted.mean(dalys, person_days),
      person_days = sum(person_days),
      time = mean(time),
      .by = c(year, age_lower, age_upper)
    ) %>%# add identifying information to output
    mutate(scen_name = scen_name_val,
           pfpr = pfpr_val,
           seas_name = seas_name_val,
           parameter_draw = run_parameters$parameter_draw,
           population = run_parameters$population,
           burnin = run_parameters$burnin)

  message('calculating prevalence')
  # Get prevalence
  prev <- raw_output %>%
    postie::get_prevalence(
      diagnostic = 'pcr',
      baseline_year = 1
    )

  prev_annual <- prev %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), mean),
      time = mean(.data$time),
      .by = c('year')
    )

  daily_output <- left_join(rates, prev)
  annual_output <- left_join(rates_annual, prev_annual)

  # Get infectivity
  message('calculating infectivity')
  infectivity <- raw_output %>%
    select(timestep,
           infectivity, infectivity_under5, infectivity_SAC, infectivity_16plus,
           starts_with('n_age')) %>%
    mutate(year = ceiling(timestep / 365),
           time = timestep) %>%
    mutate(n_age_0_36500 = 10000,
           prop_under5 = n_age_0_1825 / n_age_0_36500,
           prop_SAC = n_age_1826_5840 / n_age_0_36500,
           prop_16plus = n_age_5841_36500 / n_age_0_36500) %>%
    rowwise() %>%
    mutate(infectivity_total = infectivity_under5 + infectivity_SAC + infectivity_16plus) %>%
    ungroup() %>%
    mutate(prop_inf_under5 = infectivity_under5 / infectivity_total,
           prop_inf_SAC = infectivity_SAC / infectivity_total,
           prop_inf_16plus = infectivity_16plus / infectivity_total) %>%# add identifying information to output
    mutate(scen_name = scen_name_val,
           pfpr = pfpr_val,
           seas_name = seas_name_val,
           parameter_draw = run_parameters$parameter_draw,
           population = run_parameters$population,
           burnin = run_parameters$burnin)

  infectivity_annual <- infectivity %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), mean),
      time = mean(.data$time),
      .by = c('year')
    ) %>% select(-time, -timestep)%>%# add identifying information to output
    mutate(scen_name = scen_name_val,
           pfpr = pfpr_val,
           seas_name = seas_name_val,
           parameter_draw = run_parameters$parameter_draw,
           population = run_parameters$population,
           burnin = run_parameters$burnin)

  return(list('raw_output' = raw_output,
              'infectivity' = infectivity,
              'infectivity_annual' = infectivity_annual,
              'daily_output' = daily_output,
              'annual_output' = annual_output,
              "parameters" = params_scenario))
}


