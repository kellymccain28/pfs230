
# Defined at top level so each task is self-contained
sim_once_basic <- function(repnum, df_fakesite, par) {
  infectivity      <- par$infectivity
  offset_sd        <- par$offset_sd
  effect_size      <- par$effect_size_vals

  df <- df_fakesite %>%
    # Add infectivity value to dataset
    mutate(infectivity = infectivity) %>%
    # Add vaccine effect for TBV group
    mutate(infectivity = ifelse(arm == 'tbv', infectivity * (1-effect_size), infectivity)) %>%
    # Convert to log odds scale
    mutate(logodds = log(infectivity / (1 - infectivity)))

  # Random noise for heterogeneity between individuals
  df$offset_participant <- rnorm(nrow(df), 0, sd = offset_sd)

  # Add multiple visits per person and 60 mosquitoes per person
  d_mosq <- crossing(df, visit = c(1,4), mosquitoes = 1:60) %>%
    # Make participant-visit variable
    mutate(par_visit = paste0(ID, '_', visit)) %>%
    # Add random noise for heterogeneity between person-visits (sd is same as participant random variability, for now)
    group_by(visit, ID) %>%
    mutate(offset_participant_visit = rnorm(1, 0 , sd = offset_sd))

  d_mosq <- d_mosq %>%
    # Add offsets to the logodds and convert back to probability scale
    mutate(infectivity = plogis(logodds + offset_participant + offset_participant_visit))  # 1 / (1 + exp(-log_odds)) or inverse logit

  # Binomial draw for each row
  d_mosq$positive <- rbinom(nrow(d_mosq), size = 1, prob = d_mosq$infectivity)

  model <- glmer(positive ~ arm + (1 | ID) + (1 | par_visit),
                 data = d_mosq, family = 'binomial')
  info <- coef(summary(model))["armtbv", c("Estimate", "Pr(>|z|)")]

  data.frame(estimate = info[1],
             pval = info[2],
             sig = as.integer(info[2] < 0.05),
             correct_dir = as.integer(info[1] < 0),
             offset_sd = offset_sd,
             effect_size = effect_size,
             infectivity = infectivity,
             row.names = NULL)
}

run_power_calc_basicinf <- function(df_fakesite, nsim, pars, path_to_save, seed = 123) {
  library(tidyverse)
  library(lme4)

  cl <- parallel::makeCluster(as.integer(Sys.getenv("CCP_NUMCPUS")), outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  invisible(parallel::clusterCall(cl, ".libPaths", .libPaths()))
  invisible(parallel::clusterCall(cl, function() {
    library(tidyverse); library(lme4); TRUE
  }))
  parallel::clusterExport(cl, c("sim_once_basic", "nsim", "pars"), envir = environment())
  parallel::clusterSetRNGStream(cl, seed)   # independent, reproducible streams

  # One task per parameter set
  tasks <- lapply(pars, function(p) list(df_fakesite = df_fakesite, par = p)
  )

  # Load-balanced: idle workers pick up the next task
  results <- parallel::clusterApplyLB(cl, tasks, function(task) {
    message("running",
            " mult=", task$par$infectivity,
            " eff=", task$par$effect_size_vals,
            " sd=", task$par$offset_sd)
    bind_rows(lapply(seq_len(nsim), sim_once_basic, df_fakesite = task$df_fakesite, par = task$par))
  })

  results_df <- bind_rows(results)
  saveRDS(results_df, file.path(path_to_save, paste0('power_infranges_output_', nsim, 'sims_4visits.rds')))
  results_df
}

