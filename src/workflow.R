library(hipercow)
library(tidyverse)

# Prep hipercow:
hipercow_environment_create(sources = c("src/helper_functions.R",
                                        "src/create_site_file.R"))
# hipercow_provision(method = 'script')



# Run simulations for all sites

# List of sites to use
# data.frame with country, name of admin1 unit
countries <- c('BEN','BFA','GHA','KEN','MLI','TZA')
admin1s <- c('Atlantique','Centre-Sud','Greater Accra','Kisumu','Koulikoro','Pwani')
ur <- 'rural'
set.seed(2805)
drawnums <- sample(1:1000, 50)
parameter_draws <- c(531, 884, 125, 556,  66, 709, 876, 873, 388, 856,
                     226, 176, 272, 280, 319, 261, 432, 685, 818, 506,
                     729, 993, 296, 166, 836, 985, 936, 115, 510, 353)#,
                     # 182, 521, 793, 887, 349, 241, 765, 252, 425, 257,
                     # 270, 356,  83, 933,  90, 766, 340,  31, 324, 169)
site_df <- data.frame(country_code = countries,
                      admin_1_name = admin1s,
                      ur = 'rural')
site_df <- crossing(site_df, parameter_draw= parameter_draws)

# Only for MAP prevalence ranges
ranges <- c('lower','upper','central')
site_df_MAP <- crossing(site_df, ranges) %>%
  mutate(key = paste(country_code, admin_1_name, ur, ranges, parameter_draw, sep = '_'))
site_list_MAP <- split(site_df_MAP, seq(nrow(site_df_MAP)))

# For Jen's data calibrated EIRs
########################## DOn't add this if I want to run central/lower/upper
site_df <- site_df %>%
  mutate(ranges ='data_calibrated_weighted_MAP_testnewinfectivity',#data_calibrated
         parasit_calibration = 'weightedvar_MAP',#'weightedvar',#'weightedse',
         key = paste(country_code, admin_1_name, ur, ranges, parameter_draw, sep = '_'))
site_list_datacalib <- split(site_df, seq(nrow(site_df)))

# site_list <- site_list_MAP
site_list <- site_list_datacalib

# Locally, sequentially
source('src/create_site_file.R')
create_site_file(site_list,
                 path_to_save = paste0('outputs/jen_map_test_new_infectivity/'))


# With cluster (create_site_file.R is just a wrapper to run the analysis )
cores <- if(length(site_list) <= 32) length(site_list) else 32
t1 <- task_create_expr(expr = create_site_file(site_list,
                                               path_to_save = paste0('outputs/jen_map_test_new_infectivity_cluster/')),
                       resources = hipercow_resources(cores = cores))
task_log_show(t1)



# Combine and summarize outputs over the parameter draws
source('M:/Kelly/postdoc_JoeC/pfs230/src/combine_summarize_dfs.R')
combine_summarize_dfs(path = 'outputs/jen_map_test_new_infectivity_cluster/model_outputs/')
# combine_summarize_dfs(path = 'outputs/2026-09-08/model_outputs/')


# Plot infectivity daily, monthly, annually (for summarized runs)
path <- 'outputs/2026-09-08/model_outputs/'
outputs_processed <- readRDS(paste0(path, "processed_output_summ.rds")) # change depending on where this is saved (by date now)
# difference with below is that the outputs processed have all sites in single dfs, whereas below, they are in lists for each site (below is a precursor)
lapply(admin1s, function(a1){
  plot_infectivity_summ(a1,
                        outputs_processed,
                        time_unit = 'annual',
                        path_to_save = path)
  message('plotted annual infectivity for ', a1)

  plot_infectivity_summ(site_name,
                        outputs_processed,
                        time_unit = 'monthly',
                        path_to_save = path)
  message('plotted monthly infectivity for ', a1)

  plot_infectivity_summ(site_name,
                        outputs_processed,
                        time_unit = 'daily',
                        path_to_save = path)
  message('plotted daily infectivity for ', a1)

})

# Plot infectivity daily, monthly and annually (for the single runs)
path <- 'outputs/2026-09-07/'
outputs_processed1 <- readRDS(paste0(path, "all_processed_output.rds")) # change depending on where this is saved (by date now)

lapply(outputs_processed, function(op){
  plot_infectivity(op,
                   time_unit = 'annual',
                   path_to_save = path)
  message('plotted annual infectivity for ', op$model_input$site_name, ' ', op$model_input$target_type)

  plot_infectivity(op,
                   time_unit = 'daily',
                   path_to_save = path)
  message('plotted daily infectivity for ', op$model_input$site_name, ' ', op$model_input$target_type)

  plot_infectivity(op,
                   time_unit = 'monthly',
                   path_to_save = path)
  message('plotted monthly infectivity for ', op$model_input$site_name, ' ', op$model_input$target_type)

  message('finished ', op$model_input$site_name, ' ', op$model_input$target_type)

})
