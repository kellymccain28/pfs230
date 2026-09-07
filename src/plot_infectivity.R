
plot_infectivity <- function(processed_output,
                             time_horizon){

  key <- paste0(processed_output$country, '_', processed_output$site_name, '_', processed_output$ur)


  #### Plot the results
  infectivity_all <- purrr::map_df(results, "infectivity", .id = "parameter_draw") %>%
    rename(propinf_under5 = prop_inf_under5,
           propinf_SAC = prop_inf_SAC,
           propinf_16plus = prop_inf_16plus,
           infectivity_orig = infectivity) %>%
    pivot_longer(cols = c(prop_under5, prop_SAC, prop_16plus,
                          infectivity_under5, infectivity_SAC, infectivity_16plus,
                          propinf_under5, propinf_SAC, propinf_16plus),
                 names_to = c('.value','age_grp'),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(pfpr = factor(scales::percent(pfpr),
                         levels = scales::percent(unique(pfpr))))

  infectivity_summ <- infectivity_all %>%
    filter(year == 5) %>%
    group_by(scen_name, pfpr, seas_name, parameter_draw, age_grp) %>%
    summarise(prop_med = median(prop),
              propinf_med = median(propinf),
              infectivity_med = median(infectivity),
              prop_max = max(prop),
              propinf_max = max(propinf),
              infectivity_max = max(infectivity),
              prop_min = min(prop),
              propinf_min = min(propinf),
              infectivity_min = min(infectivity)

    )

  # proportion infectivity in each
  p1 <- ggplot(infectivity_summ) +
    geom_point(aes(x = pfpr,
                   y = propinf_med,
                   color = age_grp,
                   group = age_grp)) +
    geom_line(aes(x = pfpr,
                  y = propinf_med,
                  color = age_grp,
                  group = age_grp)) +
    facet_wrap(~seas_name) +
    labs(x = 'PfPR',
         y = 'Proportion infectivity',
         color = 'Age group',
         subtitle = key) +
    theme_classic(base_size = 12)

  # proportion of total infectivity in time horizon
  p2 <- ggplot(infectivity_summ) +
    geom_col(aes(x = pfpr,
                 y = propinf_med,
                 fill = age_grp,
                 group = age_grp),
             position = position_fill()) +
    facet_wrap(~seas_name) +
    labs(x = 'PfPR',
         y = 'Proportion infectivity',
         fill = 'Age group',
         subtitle = key) +
    theme_classic(base_size = 12)

  # absolute infectivity in time horizon
  p3 <- ggplot(infectivity_summ) +
    geom_col(aes(x = pfpr,
                 y = infectivity_med,
                 fill = age_grp,
                 group = age_grp)) +
    facet_wrap(~seas_name) +
    labs(x = 'PfPR',
         y = 'Absolute infectivity',
         fill = 'Age group',
         subtitle = key) +
    theme_classic(base_size = 12)
}
