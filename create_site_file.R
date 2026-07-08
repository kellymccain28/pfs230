## Script to run site files

# pak::pak("mrc-ide/site")
# pak::pak('mrc-ide/netz')
# devtools::install_github('mrc-ide/malariasimulation@Pfs230_2026')
create_site_file <- function(){
  library(malariasimulation)
  library(ggplot2)
  library(tidyverse)
  library(site)
  library(cowplot)
  library(netz)
  library(hipercow)

  year <- 365
  source('helper_functions.R')

  # List of sites to use
  # data.frame with country, name of admin1 unit
  countries <- c('BEN','BFA','GHA','KEN','MLI','TZA')
  admin1s <- c('Atlantique','Centre-Sud','Greater Accra','Kisumu','Koulikoro','Pwani')
  ur <- 'rural'
  site_df <- data.frame(country_code = countries,
                        admin_1_name = admin1s,
                        ur = 'rural')
  site_list <- split(site_df, seq(nrow(site_df)))

  #### Fetch site files for all sites: to be run once -- needs to be run again if site_list is updated
  # produces saved site files, individually and combination
  # site_files <- purrr::pmap(site_list, function(country_code, admin_1_name, ur){
  #   fetch_all_sites(country_code = country_code,
  #                   admin_1_name = admin_1_name,
  #                   ur = ur)
  # })
  #
  # names(site_files) <- paste0(site_list$country_code, '_', site_list$admin_1_name, '_', site_list$ur)
  # saveRDS(site_files, 'site_files/all_site_files.rds')


  # Run analysis code: ----
  cluster_cores <- Sys.getenv("CCP_NUMCPUS")


  if (cluster_cores == "") {

    message("running in serial")

    lapply(site_list,
           function(s){
             run_analysis(site = s,
                          quick_run = TRUE,
                          parameter_draw = 0)
           })
  } else {

    message(sprintf("running in parallel on %s (on the cluster)", cluster_cores))

    cl <- parallel::makeCluster(as.integer(cluster_cores),
                                outfile ="")

    invisible(parallel::clusterCall(cl, ".libPaths", .libPaths()))
    message('stop 4')
    parallel::clusterCall(cl, function() {
      message('running')
      library(malariasimulation)
      library(ggplot2)
      library(tidyverse)
      library(site)
      library(cowplot)
      library(netz)
      library(hipercow)
      library(reshape2)

      source('M:/Kelly/postdoc_JoeC/pfs230/helper_functions.R')

      TRUE
    })

    parallel::clusterExport(cl, "site_list",
                            envir = environment())

    message('stop 5')

    parallel::clusterApply(cl,
                           site_list,
                           function(s){
                             run_analysis(site = s,
                                          quick_run = TRUE,
                                          parameter_draw = 0)
                           })

    parallel::stopCluster(cl)

  }
}

# for 1 site:
# run_analysis(site = site_list[[1]],
#              quick_run = TRUE,
#              parameter_draw = 0)







# Below is the code to run outside of run_analysis function ----

# Then run gather_params for the whole list of sites above
# which should extract important bits of ach of the site files and add in ages,etc.
# model_input <- gather_params(site_files$TZA_Pwani_rural,
#                         hum_pop = 20000,
#                         quick_run = TRUE,
#                         parameter_draw = 0)

# all_model_input <- lapply(site_files,
#                      gather_params,
#                      hum_pop = 20000,
#                      quick_run = TRUE)
# names(all_model_input) <- paste0(site_list$country_code, '_', site_list$admin_1_name, '_', site_list$ur)
# saveRDS(all_model_input, 'site_files/all_model_input.rds')


# Calibration of the model for each site


# output of above used as input into a run_model function
# output <- run_model(model_input = model_input,
#                     verbose = TRUE)
# saveRDS(output, 'testoutput.rds')
#
# all_model_output <- lapply(all_model_input,
#                            run_model)
# names(all_model_output) <- paste0(site_list$country_code, '_', site_list$admin_1_name, '_', site_list$ur)
# saveRDS(all_model_output, 'all_model_output.rds')
#
#
# # process the model output -- i.e. with postie as in infectivity work from the other day
# site_list$site_key <- paste0(site_list$country_code, '_', site_list$admin_1_name, '_', site_list$ur)
# all_processed_output <- list()
# for (s in site_list$site_key){
#   message(s)
#   output_processed <- process_output(all_model_output[[s]],
#                                      all_model_input[[s]])
#
#   all_processed_output[[s]] <- output_processed
#
# }
# names(all_processed_output) <- paste0(site_list$country_code, '_', site_list$admin_1_name, '_', site_list$ur)
# saveRDS(all_processed_output, 'outputs/all_processed_output.rds')
#
#
#
# # plot the interventions over time using postie
# for (s in site_list$site_key){
#   message(s)
#
#   plot_site_files(model = all_model_output[[s]],
#                   model_input = all_model_input[[s]],
#                   site_file = site_files[[s]])
#
# }
#
#
# # plot the cases/infectivity/etc.
# for(s in site_list$site_key){
#   message(s)
#
#   plot_infectivity(all_processed_output[[s]],
#                    time_unit = 'annual')
#
#   plot_infectivity(all_processed_output[[s]],
#                    time_unit = 'daily')
# }

