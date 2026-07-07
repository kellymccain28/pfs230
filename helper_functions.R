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


# Diagnostics for site files
# function to look at the interventions over time
# Function to plot interventions over time

#' @param site_file site file for 1 site
#' @param model raw model output
#' @param model_input run parameters and parameter list
plot_site_files <- function(model,
                            model_input,
                            site_file){

  params <- model_input$param_list
  key <- paste0(site_file$sites$country, '_', site_file$sites$name_1, '_', site_file$sites$urban_rural)

  #equal spacing around colour wheel
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  #Treatment coverage / drugs used
  dr <- data.frame('t'= 1985+(params$clinical_treatment_timesteps[[1]]-1)/365,
                   'non_act' =  params$clinical_treatment_coverages[[1]],
                   'act' = params$clinical_treatment_coverages[[2]])
  drm <- reshape2::melt(dr, id.vars = 't')

  n = 2
  cols = gg_color_hue(n)
  pl1 <- ggplot(drm) +
    geom_col(aes(x = t, y = value, fill = variable),
             position = position_stack(reverse = TRUE)) +
    theme_classic() + xlim(c(1999.5,NA)) +
    labs(y = 'Treatment Coverage',
         x = 'Year',
         title = key) +
    ylim(c(0,1)) + theme(legend.position = c(0.2,0.78)) +
    scale_fill_manual(values = cols, labels = c('Non-ACT','ACT'),
                      name = 'Treatment type')
  pl1

  params$bednet_coverages
  params$bednet_rn

  #if you've run the model
  dfi <- data.frame('yr' = 1985 + model$timestep/365,
                    'itn_use' = model$n_use_net/model$n_age_0_36500[9000])
  tst <- data.frame('yr' = site_file$interventions$itn$use$year + 0.5,
                    'itn_use' = site_file$interventions$itn$use$itn_use,
                    'itn_type' = site_file$interventions$itn$implementation$net_type)
  pl2 <- ggplot(dfi) + geom_line(aes(x = yr, y = itn_use)) +
    geom_point(data = tst, aes(x = yr, y = itn_use, color = itn_type)) +
    theme_classic() +
    labs(y = 'ITN use in the population',
         x = 'Year',
         title = key) +
    scale_color_manual(values = cols, labels = c('Pyrethroid only','Pyrethroid-PBO'),
                       name = 'ITN type distributed') +
    xlim(c(1999.5,NA)) + theme(legend.position = c(0.15,0.7))#theme(legend.position = 'top')
  pl2

  ##Pyrethroid resistance
  dfr <- data.frame('yr' = site_file$vectors$pyrethroid_resistance$year + 0.5,
                    'res' = site_file$vectors$pyrethroid_resistance$pyrethroid_resistance)
  pl3 <- ggplot(dfr) + geom_point(aes(x = yr, y = res), color = 'purple') +
    theme_classic() +
    xlab('Year') + ylab('Degree of pyrethroid resistance') + ylim(c(0,1))
  pl3

  #SMC

  dfs <- data.frame('yr' = site_file$interventions$smc$implementation$year + 0.5,
                    'smc_cov' = site_file$interventions$smc$implementation$smc_cov,
                    'smc_rounds' = site_file$interventions$smc$implementation %>% ungroup() %>%
                      group_by(year) %>% mutate(n_rounds = max(round)) %>% pull(n_rounds))
  pl4 <- ggplot(dfs) + geom_point(aes(x = yr, y = smc_cov, color = factor(smc_rounds))) +
    theme_classic() +
    labs(y = 'SMC coverage',
         x = 'Year',
         title = key) +
    theme(legend.position = c(0.2,0.78)) + ylim(c(0,1)) +
    scale_color_manual(values = c('slateblue'),
                       name = 'No. of\nrounds\nper year')
  pl4

  #species & seasonality
  dsp <- data.frame('yr' = 1985 + model$timestep/365,
                    'fun'= model$total_M_funestus/model$n_age_0_36500[9000],
                    'gamb' = model$total_M_gambiae/model$n_age_0_36500[9000],
                    'arab' = model$total_M_arabiensis/model$n_age_0_36500[9000],
                    'inc' = 18*model$n_inc_clinical_1825_5840,
                    'totM' = (model$total_M_funestus+model$total_M_gambiae+model$total_M_arabiensis)/model$n_age_0_36500[9000])
  dspm <- reshape2::melt(dsp, id.vars = 'yr')
  ggplot(dspm) + geom_line(aes(x = yr, y = value, color = variable)) +
    theme_classic() + xlim(c(2010,2012))

  seas <- function(model){
    mn <- min(model$timestep)
    mx <- max(model$timestep)
    store <- rep(0,365)
    for(i in mn:mx){
      sr <- model$timestep[i] %% 365
      store[sr] <- store[sr] + model$n_inc_clinical_0_36500[i]/model$n_age_0_36500[i]
    }
    bray <- cumsum(store)
    return(bray/bray[365])
  }

  fg <- seas(model = model)
  dfg <- data.frame('t' = seq(1,365,1), fg = fg)
  pl5 <- ggplot(dfg) + geom_hline(yintercept = 100*c(0.1,0.9), color = 'red', alpha = .3) +
    geom_line(aes(x=t,y=100*fg)) + theme_classic() +
    scale_x_continuous(breaks = c(0,91,182,274),
                       labels = c('Jan','April','July','Oct')) +
    labs(y = 'Averaged cumulative incidence (%)',
         title = paste0(key, ': Seasonality of malaria incidence'),
         x = 'Month')
  pl5

  # plaux <- cowplot::plot_grid(pl1,pl3,pl4,pl5, nrow = 2)
  # cowplot::plot_grid(pl2,plaux, nrow = 2, rel_heights = c(0.6,1))
  # ggsave('pfVIMT_site_files/BE.pdf', height = 8.7, width = 9)

  # Open the PDF device and specify the file path and name
  pdf(file = paste0("site_files/site_file_plots/", key, '.pdf'))

  # Generate plots
  print(pl1)
  print(pl2)
  print(pl3)
  print(pl4)
  print(pl5)

  # Close the PDF device to finalize the file
  dev.off()
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
    term_yr<- 2026

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

# Calibration of the model to site-specific information -- to check later
pr_match <- function(x, y){

  data <- baseline_parameters[x,]
  params <- unlist(data$params, recursive = FALSE)
  params$timesteps <- data$sim_length + data$warmup

  # defining target as pfpr value in last 2 years of simulation
  target <- rep(data$pfpr, 2)

  set.seed(1234)
  out <- cali::calibrate(parameters = params,
                         target = target,
                         summary_function = annual_pfpr_summary,
                         tolerance = 0.001,
                         low = 0.1,
                         high = 1500)

  # store init_EIR results as an .rds file to be read in later
  PR <- data.frame(scenarioID = x,  drawID = y)
  PR$starting_EIR <- out
  PR$ID <- data$ID

  print(paste0('Finished scenario ',x))
  saveRDS(PR, paste0('PrEIR/PRmatch_draws_', data$ID, '.rds'))
}


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
process_output <- function(model, model_input){
  key <- paste0(model_input$country, '_', model_input$site_name, '_', model_input$ur)

  # Drop burnin
  raw_output <- postie::drop_burnin(model, burnin = model_input$burnin * 365)

  message('calculating rates')

  rates <- postie::get_rates(
    raw_output,
    baseline_year = model_input$param_list$start_year + model_input$burnin
  ) %>%# add identifying information to output
    mutate(country = model_input$country,
           ur = model_input$ur,
           site_name = model_input$site_name,
           parameter_draw = model_input$parameter_draw,
           population = model_input$population,
           burnin = model_input$burnin)

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
    mutate(country = model_input$country,
           ur = model_input$ur,
           site_name = model_input$site_name,
           parameter_draw = model_input$parameter_draw,
           population = model_input$population,
           burnin = model_input$burnin)

  message('calculating prevalence')
  # Get prevalence
  prev <- raw_output %>%
    postie::get_prevalence(
      diagnostic = 'pcr',
      baseline_year = model_input$param_list$start_year + model_input$burnin
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
    mutate(year = floor(timestep / 365) + model_input$param_list$start_year + model_input$burnin,
           time = timestep) %>%
    mutate(prop_under5 = n_age_0_1825 / n_age_0_36500,
           prop_SAC = n_age_1825_5840 / n_age_0_36500,
           prop_16plus = n_age_5840_36500 / n_age_0_36500) %>%
    rowwise() %>%
    # get total summed infectivity
    mutate(infectivity_sum_total = infectivity_under5 + infectivity_SAC + infectivity_16plus) %>%
    ungroup() %>%

    # Proportion of summed infectivity of total summed infectivity
    mutate(prop_sum_inf_under5 = infectivity_under5 / infectivity_sum_total,
           prop_sum_inf_SAC = infectivity_SAC / infectivity_sum_total,
           prop_sum_inf_16plus = infectivity_16plus / infectivity_sum_total) %>%

    # get mean per-person infectivity for each age group
    mutate(mean_inf_under5 = infectivity_under5 / n_age_0_1825,
           mean_inf_SAC = infectivity_SAC / n_age_1825_5840,
           mean_inf_16plus = infectivity_16plus / n_age_5840_36500) %>%

    rowwise() %>%
    # sum of means above should be approx equal to the originaly infectvity output
    mutate(sum_mean_infectivity = mean_inf_under5 + mean_inf_SAC + mean_inf_16plus) %>%
    ungroup() %>%

    # get proportion of mean infectivity by age group
    mutate(prop_mean_inf_under5 = mean_inf_under5 / sum_mean_infectivity,
           prop_mean_inf_SAC = mean_inf_SAC / sum_mean_infectivity,
           prop_mean_inf_16plus = mean_inf_16plus / sum_mean_infectivity) %>%

    # add identifying information to output
    mutate(country = model_input$country,
           ur = model_input$ur,
           site_name = model_input$site_name,
           parameter_draw = model_input$parameter_draw,
           population = model_input$population,
           burnin = model_input$burnin)

  infectivity_annual <- infectivity %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), mean),
      time = mean(.data$time),
      .by = c('year')
    ) %>% select(-time, -timestep)%>%# add identifying information to output
    mutate(country = model_input$country,
           ur = model_input$ur,
           site_name = model_input$site_name,
           parameter_draw = model_input$parameter_draw,
           population = model_input$population,
           burnin = model_input$burnin)

  processed_out <- list('raw_output' = raw_output,
                        'infectivity' = infectivity,
                        'infectivity_annual' = infectivity_annual,
                        'daily_epi_output' = daily_output,
                        'annual_epi_output' = annual_output,
                        "model_input" = model_input)

  if(!dir.exists(paste0("outputs/", key, "/"))){
    dir.create(paste0("outputs/", key, "/"))
  }

  saveRDS(processed_out, paste0('outputs/', key, '/processed_out_', key, '.rds'))

  return(processed_out)
}


# Plot infectivity ----
#' @param processed_output processed output for a single site
plot_infectivity <- function(processed_output,
                             time_unit){

  key <- paste0(processed_output$model_input$country, '_', processed_output$model_input$site_name, '_', processed_output$model_input$ur)


  if(time_unit == 'annual'){
    inf <- processed_output$infectivity_annual
    inf$time <- inf$year
  } else if(time_unit == 'daily'){
    inf <- processed_output$infectivity
    inf$time <- inf$timestep / 365
  }

  # ggplot(inf) +
  #   geom_point(aes(x = time, y = infectivity))

  inf_long <- inf %>%
    select(time, infectivity_under5, infectivity_SAC, infectivity_16plus,
           mean_inf_under5, mean_inf_SAC, mean_inf_16plus,
           prop_mean_inf_under5, prop_mean_inf_SAC, prop_mean_inf_16plus) %>%
    pivot_longer(
      cols = -time,
      names_to = c(".value", "age_group"),
      names_pattern = "(infectivity|mean_inf|prop_mean_inf)_(under5|SAC|16plus)"
    )

  p1 <- ggplot(inf_long) +
    geom_line(aes(x = time, y = infectivity, color = age_group)) +
    labs(y = 'Infectivity sum by age group',
         x = 'Year',
         title = key) +
    theme_classic(base_size = 12)

  p2 <- ggplot(inf_long) +
    geom_line(aes(x = time, y = mean_inf, color = age_group)) +
    labs(y = 'Mean infectivity per person',
         x = 'Year',
         title = key) +
    theme_classic(base_size = 12)

  p3 <- ggplot(inf_long) +
    geom_line(aes(x = time, y = prop_mean_inf, color = age_group)) +
    labs(y = 'Proportion of mean infectivity per person',
         x = 'Year',
         title = key) +
    theme_classic(base_size = 12)

  # At last timestep or last year
  infectivity_summ <- inf %>%
      filter(time == max(inf$time)) %>%
      select(time, infectivity_under5, infectivity_SAC, infectivity_16plus,
             mean_inf_under5, mean_inf_SAC, mean_inf_16plus,
             prop_mean_inf_under5, prop_mean_inf_SAC, prop_mean_inf_16plus) %>%
      pivot_longer(
        cols = -time,
        names_to = c(".value", "age_group"),
        names_pattern = "(infectivity|mean_inf|prop_mean_inf)_(under5|SAC|16plus)"
      )

  p4 <- ggplot(infectivity_summ) +
    geom_col(aes(x = age_group, y = prop_mean_inf), fill = 'darkred') +
    labs(y = 'Proportion of mean infectivity per person',
         x = 'Year',
         title = paste0(key,' ', ifelse(time_unit == 'annual', ' in last timestep', 'in last year')))+
    theme_classic(base_size = 12)

    # Open the PDF device and specify the file path and name
  if(!dir.exists(paste0("outputs/", key, "/"))){
    dir.create(paste0("outputs/", key, "/"))
  }
  pdf(file = paste0("outputs/", key, "/infectivity", key, '_', time_unit, '.pdf'))

  # Generate plots
  print(p1)
  print(p2)
  print(p3)
  print(p4)

  # Close the PDF device to finalize the file
  dev.off()

}
