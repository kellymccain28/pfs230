# With cluster (run_power_calc.R is just a wrapper to run the analysis )
# site_specific must be a list of trial pops per site
# nsim is the number of repetions to run
run_power_calc <- function(site_specific, nsim, pars, path_to_save){
  library(tidyverse)
  library(lme4)
  message('loaded packages')
  # Run power analysis code: ----
  cluster_cores <- Sys.getenv("CCP_NUMCPUS")
  cl <- parallel::makeCluster(as.integer(cluster_cores),
                              outfile ="")
  message('created cluster')
  invisible(parallel::clusterCall(cl, ".libPaths", .libPaths()))

  parallel::clusterCall(cl, function() {
    message('running')
    library(tidyverse)
    library(lme4)

    source('M:/Kelly/postdoc_JoeC/pfs230/power_calculation/run_power_calc.R')

    TRUE
  })

  parallel::clusterExport(cl, c("site_specific","nsim",
                                "pars"),
                          envir = environment())
  message('exported items')
  results <- parallel::clusterApply(cl,
                                    site_specific,
                                    function(df_site, nsims = nsim){
                                      message("running site ", df_site$country[1], " ", nsims, " times")

                                      sim_once <- function(repnum, par){ # rep number and row with things i've varied
                                        infectivity_mult <- par$infectivity_mult
                                        offset_sd <- par$offset_sd
                                        effect_size <- par$effect_size_vals

                                        df <- df_site %>%
                                          # Apply multiplier to infectivity starting point
                                          mutate(mean_inf_median = infectivity_mult * mean_inf_median) %>%
                                          # Add vaccine effect for TBV group
                                          mutate(mean_inf_median = ifelse(arm == 'tbv', mean_inf_median * (1-effect_size), mean_inf_median)) %>%
                                          # Convert to log odds scale
                                          mutate(logodds = log(mean_inf_median / (1 - mean_inf_median))) %>%
                                          # Random noise for heterogeneity between individuals
                                          rowwise() %>% mutate(offset_participant = rnorm(1, 0, sd = offset_sd)) %>% ungroup()

                                        # Add multiple visits per person and 60 mosquitoes per person
                                        d_mosq <- crossing(df, visit = c(1,2), mosquitoes = 1:60) %>%
                                        # Add random noise for heterogeneity between person-visits
                                          group_by(visit, ID) %>%
                                          mutate(offset_participant_visit = rnorm(1, 0 , sd = offset_sd)) %>% ungroup() %>%
                                          # Add offset to the logodds
                                          mutate(logodds_new = logodds + offset_participant + offset_participant_visit) %>%
                                          # Convert back to probability scale
                                          mutate(infectivity = plogis(logodds_new))  # 1 / (1 + exp(-log_odds)) or inverse logit

                                        # Binomial draw for each row
                                        d_mosq$positive <- rbinom(nrow(d_mosq), size = 1, prob = d_mosq$infectivity)
                                        # table(d_mosq$positive, d_mosq$arm)


                                        # GLM - with random effect for participant id and for visit number
                                        model <- glmer(positive ~ arm + (1|ID) + (1|visit),
                                                       data = d_mosq,
                                                       family = 'binomial')
                                        info <- coef(summary(model))["armtbv", c("Estimate","Pr(>|z|)")]

                                        o <- data.frame(estimate = info[1],
                                                        pval = info[2],
                                                        sig = ifelse(info[2] < 0.05, 1, 0),
                                                        correct_dir = ifelse(info[1] < 0, 1, 0),
                                                        country = df_site$country[1],
                                                        rep = repnum,
                                                        offset_sd = offset_sd,
                                                        effect_size = effect_size,
                                                        infectivity_mult = infectivity_mult)
                                        rownames(o) <- NULL

                                        return(o)
                                      }

                                      bind_rows(lapply(pars, function(p) {
                                        bind_rows(lapply(1:nsims, sim_once, p))
                                      } ))

                                    })

  results_df <- bind_rows(results)

  saveRDS(results_df, paste0(path_to_save, '/power_output_', nsim, 'sims.rds'))

  parallel::stopCluster(cl)

}

