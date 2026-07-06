## Script to run site files

# pak::pak("mrc-ide/site")
# pak::pak('mrc-ide/netz')
# devtools::install_github('mrc-ide/malariasimulation@Pfs230_2026')
library(malariasimulation)
library(ggplot2)
library(tidyverse)
library(site)
library(cowplot)
library(netz)


year <- 365
source('helper_functions.R')

# List of sites to use
# data.frame with country, name of admin1 unit
countries <- c('BEN','BFA','GHA','KEN','MLI','TZA')
admin1s <- c('Atlantique','Centre-Sud','Greater Accra','Kisumu','Koulikoro','Pwani')
ur <- 'rural'
site_list <- data.frame(country_code = countries,
                        admin_1_name = admin1s,
                        ur = 'rural')

#### Fetch site files for all sites: to be run once
# produces saved site files
site_files <- purrr::pmap(site_list, function(country_code, admin_1_name, ur){
  fetch_all_sites(country_code = country_code,
                  admin_1_name = admin_1_name,
                  ur = ur)
})

names(site_files) <- paste0(site_list$country_code, '_', site_list$admin_1_name, '_', site_list$ur)
saveRDS(site_files, 'site_files/all_site_files.rds')

# Then run gather_params for the whole list of sites above
# which should extract important bits of ach of the site files and add in ages,etc.

params <- gather_params(site_files$TZA_Pwani_rural,
                        hum_pop = 20000,
                        quick_run = TRUE,
                        parameter_draw = 0)

all_params <- lapply(site_files,
                     gather_params,
                     hum_pop = 20000,
                     quick_run = TRUE)

# Calibration of the model for each site


# output of above used as input into a run_model function
output <- run_model(model_input = params,
          verbose = TRUE)
saveRDS(output, 'testoutput.rds')

# process the model output -- i.e. with postie as in infectivity work from the other day
output_processed <- process_output(output)

# plot the interventions over time using postie

# plot the cases/infectivity/etc.


