# Power calculation
library(lme4)

# Define fixed values
total_pop <- 4000
sites <- data.frame(country_code = c('BFA','BEN','KEN','MLI','TZA','GHA'),
                    country = c('Burkina Faso','Benin','Kenya','Mali','Tanzania','Ghana'))
n_sites <- nrow(sites)
n_per_site <- floor(total_pop / n_sites)
effect_size <- 0.75 # in terms of efficacy, so 0.75 means 75% efficacy or x * (1-0.75)
offset_sd <- 0.3

# Define dataset of trial site populations
# p_child <- 0.20 # proportion of pop aged 9-17 (because this pop is enriched, can change value as needed)
# sample_ages <- function(n, p_child = 0.30) {
#   n_child <- round(n * p_child)
#   n_adult <- n - n_child
#   c(sample(9:17,  n_child, replace = TRUE),   # enriched pediatric group
#     sample(c(5:8, 18:99), n_adult, replace = TRUE))  # everyone else 5-99
# }

hpop <- data.frame()
for (i in seq_along(sites$country_code)) {
  # Assign arms in a 3:1 vaccine:control ratio.
  n_vaccine <- round(n_per_site * 3 / 4)
  n_control <- n_per_site - n_vaccine
  arms      <- sample(rep(c("tbv","control"), c(n_vaccine, n_control)))

  site_df <- data.frame(
    country_code    = sites$country_code[i],
    country = sites$country[i],
    arm     = arms,
    # age     = sample_ages(n, p_child),
    stringsAsFactors = FALSE
  )
  site_df$ID <- sprintf("%s%s_%04d", site_df$country_code, site_df$arm, seq_len(n_per_site))

  hpop <- rbind(hpop, site_df)
}


# Summarize infectivity model outputs - to be used as probabilities of infectivity in bin draws
path <- 'outputs/weighted_calibrationvariance_MAP/'
infect_annual_summ <- readRDS(paste0(path, 'model_outputs/summarized_infectivity_annual.rds'))

# Get the annual infectivity data as predicted by malariasimulation
infectivity <- infect_annual_summ %>%
  filter(year == 2026) %>%
  select(country, site_name, ur, target_type, year,
         # starts_with('infectivity'),starts_with('mean_inf'), starts_with('prop_sum'),
         contains('5to17')
  ) %>%
  pivot_longer(
    cols = -c('country', 'site_name', 'ur', 'target_type', 'year'),
    names_to = c("metric", "age_group", "estimate"),
    names_pattern = "^(infectivity|prop_sum_inf|mean_inf)(?:_(under5|SAC|16plus|youngSAC|oldSAC|5to17|18plus))?_(median|lower|upper)$",
    values_to = 'value'
  ) %>%
  pivot_wider(
    names_from = c(metric, estimate),
    values_from = value,
    names_glue = "{metric}_{estimate}"
  ) %>%
  select(country, site_name, mean_inf_median)

# Add modelled infectivity to dataset (control arm)
d <- left_join(hpop, infectivity)

# Run the power calculation
site_specific <- split(d, d$country)

# hipercow_environment_create(name = 'power_env',
#                             sources = c("power_calculation/run_power_calc.R"))

# Below is to run singel effect size/offset_sd
nreps = 200
pars <- crossing(infectivity_mult = seq(0.8, 1.2, 0.2),
                 effect_size_vals = seq(0.15, 0.75, 0.2),
                 offset_sd = 0.3)
pars <- split(pars, seq(nrow(pars)))
cores <- if(length(site_specific) <= 32) length(site_specific) else 32
t200 <- task_create_expr(expr = run_power_calc(site_specific,
                                             nsim = nreps,
                                             pars,
                                             path_to_save = paste0('M:/Kelly/postdoc_JoeC/pfs230/power_calculation/outputs')),
                       environment = 'power_env',
                       resources = hipercow_resources(cores = cores))
task_log_show(t1)#20 sims
task_log_show(t200)#200 sims


# Make task bundle with mutliple parameter sets
# pars <- crossing(infectivity_mult = seq(0.8, 1.2, 0.2),
#                  effect_size_vals = seq(0.25, 0.75, 0.25),
#                  offset_sd = 0.3)
# pars <- split(pars, seq(nrow(pars)))
#
# create_power_bundle <- function(x){
#   bundle1 <- task_create_bulk_expr(
#     run_power_calc(site_specific,
#                    nsim = nreps,
#                    effect_size = effect_size_vals,
#                    offset_sd = offset_sd,
#                    infectivity_mult = infectivity_mult,
#                    path_to_save = paste0('M:/Kelly/postdoc_JoeC/pfs230/power_calculation/outputs')
#     ),
#     pars[x,],
#     environment = 'power_env',
#     resources = hipercow_resources(cores = cores),
#     bundle_name = 'power_calculations')
#
#   return(bundle1)
# }
# bundle1 <- create_power_bundle(1:nrow(pars))
# table(hipercow_bundle_status(bundle1))


# Plotting the power calculations
library(ggplot2)
out <- readRDS('power_calculation/outputs/power_output_20sims.rds')
out_summ <- bind_rows(out) %>%
  mutate(sig_correctdir = ifelse(sig == 1 & correct_dir == 1, 1, 0)) %>%
  group_by(country, offset_sd, effect_size, infectivity_mult) %>%
  summarize(n_sig_correctdir = sum(sig_correctdir),
            n_per_country = n(),
            power = n_sig_correctdir / n_per_country) %>%
  rowwise() %>%
  mutate(power_lower = binom.test(n_sig_correctdir, n_per_country)$conf.int[1],
         power_upper = binom.test(n_sig_correctdir, n_per_country)$conf.int[2])

ggplot(out_summ %>% mutate(infectivity_mult = paste0('Infectivity multiplier: ', infectivity_mult))) +
  geom_line(aes(x = effect_size, y = power, color = country)) +
  geom_ribbon(aes(x = effect_size, ymin = power_lower, ymax = power_upper, fill = country),
              alpha = 0.2) +
  geom_point(aes(x = effect_size, y = power, color = country)) +
  geom_hline(aes(yintercept = 0.8), linetype = 2) +
  facet_grid(rows = vars(country),
             cols = vars(infectivity_mult)) +
  labs(x = 'Effect size',
       y = 'Power') +
  theme_bw() +
  theme(legend.position = 'none')
