library(hipercow)

#Prep hipercow:
hipercow_environment_create(sources = c("helper_functions.R",
                                        "create_site_file.R"))
hipercow_provision(method = 'script')

t1 <- task_create_expr(expr = create_site_file(),
                 resources = hipercow_resources(cores = 6))
task_log_show(t1)
