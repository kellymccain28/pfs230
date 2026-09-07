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
site_df <- data.frame(country_code = countries,
                      admin_1_name = admin1s,
                      ur = 'rural')

# Only for MAP prevalence ranges
ranges <- c('lower','upper','central')
site_df_MAP <- crossing(site_df, ranges) %>%
  mutate(key = paste(country_code, admin_1_name, ur, ranges, sep = '_'))
site_list_MAP <- split(site_df_MAP, seq(nrow(site_df_MAP)))

# For Jen's data calibrated EIRs
########################## DOn't add this if I want to run central/lower/upper
site_df$parasit_calibration = TRUE
site_df$ranges <- 'data_calibrated'
site_list_datacalib <- split(site_df, seq(nrow(site_df)))

site_list <- site_list_MAP
site_list <- site_list_datacalib

# Locally, sequentially
source('src/create_site_file.R')
create_site_file(site_list)


# With cluster (create_site_file.R is just a wrapper to run the analysis )
cores <- if(length(site_list) <= 32) length(site_list) else 32
t1 <- task_create_expr(expr = create_site_file(site_list),
                       resources = hipercow_resources(cores = cores))
task_log_show(t1)




# Plot infectivity daily and annually
outputs_processed <- readRDS("M:/Kelly/postdoc_JoeC/pfs230/outputs/all_processed_output.rds")
names(outputs_processed) <- site_df$key

lapply(outputs_processed, function(op){
  plot_infectivity(op,
                   time_unit = 'annual')
  message('plotted annual infectivity for ', op$model_input$site_name, ' ', op$model_input$target_type)

  plot_infectivity(op,
                   time_unit = 'daily')
  message('plotted daily infectivity for ', op$model_input$site_name, ' ', op$model_input$target_type)

  message('finished ', op$model_input$site_name, ' ', op$model_input$target_type)

})
