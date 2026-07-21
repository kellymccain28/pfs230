library(hipercow)

#Prep hipercow:
hipercow_environment_create(sources = c("helper_functions.R",
                                        "create_site_file.R"))
hipercow_provision(method = 'script')




# Run simulations for all sites

# List of sites to use
# data.frame with country, name of admin1 unit
countries <- c('BEN','BFA','GHA','KEN','MLI','TZA')
admin1s <- c('Atlantique','Centre-Sud','Greater Accra','Kisumu','Koulikoro','Pwani')
ur <- 'rural'
site_df <- data.frame(country_code = countries,
                      admin_1_name = admin1s,
                      ur = 'rural')
ranges <- c('lower','upper','central')
site_df <- crossing(site_df, ranges) %>%
  mutate(key = paste(country_code, admin_1_name, ur, ranges, sep = '_'))
site_list <- split(site_df, seq(nrow(site_df)))

create_site_file(site_list)

cores <- if(length(site_list) <=32) length(site_list) else 32
t1 <- task_create_expr(expr = create_site_file(site_list),
                 resources = hipercow_resources(cores = cores))
task_log_show(t1)
