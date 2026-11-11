# Script to make plots for manuscript
# uses dataset processed in malsim_tbaplus.R script
#

# =============================================================================
# Publication figures: TBV oocyst -> infection models and epidemiological impact
#
#   Fig 1  Model figure       (a) oocyst distribution +/- TBV, (b) oocysts -> P(inf),
#                             (c) TRA -> TBA+
#   Fig 2  Epi figure         (a) cumulative incidence, (b) monthly cases averted,
#                             (c) cumulative cases averted, (d) PfPR2-10, (e) infectivity
#   Fig 3  Efficiency figure  (a) monthly and (b) cumulative cases averted per 1,000 doses
#   Table  End-of-horizon summary (CSV)
#
# Figs 2 and 3 are faceted: columns = who is vaccinated (tbv_ages), rows = age group
# whose cases are counted. y-axes are shared across columns within a row so the two
# vaccination strategies can be compared directly.
#
# Inputs : summarised results .rds from the scenario script (sum_all)
# Outputs: PNG (600 dpi) + PDF (vector) for each figure in `fig_dir`
# Needs  : dplyr, tidyr, ggplot2, patchwork, scales
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---- 0. Configuration -------------------------------------------------------
input_rds <- 'stan_fitting/outputs/scenarios/summarized_allresults_pop100k_20261106_0draws_3mtbv_booster_newmodels.rds'
fig_dir   <- 'stan_fitting/outputs/figures'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

PNG_DPI      <- 600
FIG_TITLE    <- NULL      # e.g. 'Highly seasonal transmission, 80% coverage'
SHOW_CRI     <- FALSE     # TRUE once you run >1 parameter draw: adds 95% CrI ribbons
INSET_DIGITS <- 1         # decimal places (in %) for the Fig 1a inset

# Keep only scenarios whose id contains this text (so a re-run with more seasonalities /
# coverages does not mix them). Set to NULL to keep everything.
SCENARIO_MATCH <- 'highly_seasonal_cov0.8'

# tbv_ages value -> legend label (linetype)
TBV_LABELS <- c('5+'  = 'All ages 5+ vaccinated',
                'SAC' = 'School-aged children (5-15) vaccinated')
TBV_LINETYPES <- c('5+' = 'solid', 'SAC' = 'dashed')

AGE_GROUPS <- c('0-5', '0-200')        # age groups whose cases are plotted (facets)
AGE_LABELS <- c('Under 5y', 'All ages')

# ---- 1. Shared style --------------------------------------------------------
model_ids    <- c('m0', 'm1', 'm2')
lvl          <- c('none', model_ids)
model_labels <- c(none = 'No TBV',
                  m0   = 'm0',
                  m1   = 'm1',
                  m2   = 'm2')          # edit to match the manuscript wording
# Okabe-Ito (colour-blind safe); no-TBV in near-black
col_models   <- c(none = 'grey40', m0 = '#0072B2', m1 = '#D55E00', m2 = '#009E73')

scale_col_model <- function(breaks = lvl) {
  scale_colour_manual(values = col_models, breaks = breaks, limits = breaks,
                      labels = unname(model_labels[breaks]), name = NULL)
}
scale_fill_epi <- scale_fill_manual(
  values = c(col_models, 'Range across models' = 'grey55'),
  breaks = 'Range across models', name = NULL)
scale_lt_tbv <- scale_linetype_manual(
  values = TBV_LINETYPES, breaks = names(TBV_LABELS),
  labels = unname(TBV_LABELS), name = NULL)

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.border     = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = 'grey92', linewidth = 0.25),
      axis.line        = element_line(colour = 'grey20', linewidth = 0.35),
      axis.ticks       = element_line(colour = 'grey20', linewidth = 0.35),
      axis.text        = element_text(colour = 'grey15'),
      axis.title       = element_text(colour = 'grey10'),
      strip.background = element_blank(),
      strip.text       = element_text(face = 'bold', size = base_size, hjust = 0),
      panel.spacing    = unit(0.6, 'lines'),
      legend.position  = 'bottom',
      legend.key       = element_blank(),
      legend.background = element_blank(),
      legend.key.size  = unit(0.4, 'cm'),
      legend.key.width = unit(0.9, 'cm'),
      legend.text      = element_text(size = base_size),
      plot.tag         = element_text(face = 'bold', size = base_size + 4),
      plot.title       = element_text(face = 'bold', size = base_size + 2),
      plot.margin      = margin(4, 6, 4, 4)
    )
}
theme_set(theme_pub())

save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, '.png')), p, width = width, height = height,
         units = 'in', dpi = PNG_DPI, bg = 'white')
  pdf_dev <- if (isTRUE(capabilities('cairo'))) grDevices::cairo_pdf else 'pdf'
  ggsave(file.path(fig_dir, paste0(name, '.pdf')), p, width = width, height = height,
         units = 'in', device = pdf_dev)
}

# =============================================================================
# FIGURE 1: MODEL FIGURE
# =============================================================================
# The Bompard NB is so over-dispersed that ~99.998% of mosquitoes have 0 oocysts, so the
# distribution is shown for oocyst-POSITIVE mosquitoes. A TBV with TRA t thins oocysts
# binomially, which keeps the NB form: O_tbv ~ NB(mean = mu*(1-t), size = k).
# Using the analytic PMF avoids the Monte-Carlo noise of sampling 2e7 mosquitoes.
j_max  <- 3000
j_vals <- 0:j_max
j_pos  <- 1:j_max

# P(O_tbv = j | O > 0) after a TBV with TRA t (j = 0: was positive, now fully blocked)
dist_given_tra <- function(tra, mu = mALL, k = rrALL, j = j_vals) {
  mu_t  <- mu * (1 - tra)
  a     <- k * log(k / (k + mu_t))     # log P(O_tbv = 0)
  b     <- k * log(k / (k + mu))       # log P(O = 0)
  ppos  <- -expm1(b)                   # P(O > 0)
  pj    <- dnbinom(j, size = k, mu = mu_t)
  pj[j == 0] <- exp(b) * expm1(a - b)  # P(O_tbv = 0) - P(O = 0), numerically stable
  pj / ppos
}

w <- dist_given_tra(0)[-1]             # P(O = j | O > 0), j = 1..j_max
stopifnot(abs(sum(w) - 1) < 1e-6)

# Scale m1/m2 so every model has the same mean P(inf) as m0 over the oocyst distribution
m1_const <- sum(w * m0_fun(j_pos)) / sum(w * m1_raw(j_pos))
m2_const <- sum(w * m0_fun(j_pos)) / sum(w * m2_raw(j_pos))
m1_c <- function(x) m1_raw(x) * m1_const
m2_c <- function(x) m2_raw(x) * m2_const
model_funs <- list(m0 = m0_fun, m1 = m1_c, m2 = m2_c)
if (any(m1_c(j_pos) > 1) || any(m2_c(j_pos) > 1)) warning('P(inf) exceeds 1 for m1 and/or m2')

# ---- 1a. Model definitions (copied from the fitting script) ------------------
mALL    <- 0.000157     # Bompard NB mean (oocysts per mosquito, zeros included)
rrALL   <- 0.00000495   # Bompard NB dispersion (size)
pinf_b0 <- 0.590076

f2 <- function(x, beta = 10263.54, alpha = 29.95) beta * x / (1 + x / alpha)
g3 <- function(x, a = 49.61, b = 4980.93, c = 10143.51) a + b * (x^2 / (x^2 + c^2))
g6 <- function(x, floor_p = 0, mult = 0.05, scale = 34289.7) x * (floor_p + mult * 1 / (1 + x / scale))
h  <- function(x, baselineprob = 0.01, halfsat = 1000, steep = 1.5) {
  baselineprob + (pinf_b0 - baselineprob) * (x^steep / (x^steep + halfsat^steep))
}

m0_fun <- function(x, b0 = pinf_b0) ifelse(x > 0, b0, 0)
m1_raw <- function(x) ifelse(x == 0, 0, h(g3(f2(x, beta = 15000), a = 600, b = 1000, c = 20000),
                                          steep = 1.8, baselineprob = 0.05, halfsat = 3000))
m2_raw <- function(x) ifelse(x == 0, 0, h(g6(f2(x, beta = 35000), mult = 0.1, scale = 35000),
                                          steep = 7, baselineprob = 0.01, halfsat = 2200))

# ---- Figure 1 ------------------------------

# ---- Panel a: oocyst distribution from Joe ---------------------
df_NB <- readRDS('stan_fitting/outputs/pmf_tbv.rds') %>%
  rename(pmf_tra0 = pmf, cmf_tra0 = cmf) %>%   # base cols = the no-tra (notbv) case
  pivot_longer(
    cols = -o,
    names_to = c(".value", "tra"),
    names_sep = "_tra",
    names_transform = list(tra = as.integer)
  ) %>%
  mutate(
    tra = factor(tra, levels = c("0", "30", "60", "90"),
                 labels = c("No TBV", "TRA = 30%", "TRA = 60%", "TRA = 90%")))

p_a <- ggplot(df_NB) +
  geom_vline(aes(xintercept = 0), color = 'grey80', linetype = 2) +
  geom_line(aes(x = o, y = cmf, color = tra), linewidth = 0.7) +
  geom_point(aes(x = o, y = cmf, color = tra),
             shape = 21, fill = 'white', size = 1.3, stroke = 0.5) +
  theme_classic() +
  xlim(c(0,30)) +
  xlab('Oocyst count') + ylab('Cumulative mass function') +
  scale_color_manual(values = c('black','slateblue','orange','orchid'),
                     name = '', labels = c('No TBV','TRA = 30%','TRA = 60%','TRA = 90%')) +
  annotate('text', x = 14, y = df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 30%',]$cmf, colour = 'slateblue', size = 2.5,
           label = paste0(round(100*df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 30%',]$cmf),'% of transmission\nevents blocked')) +
  annotate('text', x = 14, y = df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 60%',]$cmf, colour = 'orange', size = 2.5,
           label = paste0(round(100*df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 60%',]$cmf),'% of transmission\nevents blocked')) +
  annotate('text', x = 14, y = df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 90%',]$cmf, colour = 'orchid', size = 2.5,
           label = paste0(round(100*df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 90%',]$cmf),'% of transmission\nevents blocked')) +
  annotate('segment', xend = 0.3, x = 8, y = df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 30%',]$cmf,
           'color' = 'slateblue', arrow = arrow(length = unit(0.2, "cm"), type = 'closed')) +
  annotate('segment', xend = 0.3, x = 8, y = df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 60%',]$cmf,
           'color' = 'orange', arrow = arrow(length = unit(0.2, "cm"), type = 'closed')) +
  annotate('segment', xend = 0.3, x = 8, y = df_NB[df_NB$o == 0 & df_NB$tra=='TRA = 90%',]$cmf,
           'color' = 'orchid', arrow = arrow(length = unit(0.2, "cm"), type = 'closed')) +
  theme(legend.position = c(0.85,0.4))


# ---- Panel b: oocyst count -> P(infection) -----------------------------------
x_grid  <- 1:30
pinf_df <- tibble(oocysts = x_grid,
                  m0 = m0_fun(x_grid),
                  m1 = m1_c(x_grid),
                  m2 = m2_c(x_grid)) %>%
  pivot_longer(-oocysts, names_to = 'model', values_to = 'pinf') %>%
  mutate(model = factor(model, levels = lvl))

cdf <- cumsum(w)
q05 <- j_pos[which(cdf >= 0.05)[1]]
q95 <- j_pos[which(cdf >= 0.95)[1]]

p_b <- ggplot(pinf_df, aes(x = oocysts, y = pinf, colour = model)) +
  # annotate('rect', xmin = q05, xmax = q95, ymin = -Inf, ymax = Inf, fill = 'grey85', alpha = 0.5) +
  # annotate('text', x = sqrt(q05 * q95), y = Inf, vjust = 1.3, size = 2.1, colour = 'grey30',
  #          lineheight = 0.9, label = 'central 90% of\noocyst-positive\nmosquitoes') +
  geom_point(shape = 21, fill = 'white', size = 1.3, stroke = 0.5) +
  geom_line(linewidth = 0.7) +
  scale_col_model(model_ids) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
  labs(x = 'Oocyst count (log scale)', y = 'P(infection | oocyst count)') +
  theme(legend.position = 'none')

# ---- Panel c: TRA -> TBA+ ------------------------------------------------------
# Mechanistic TBA+ for each model (exact expectation over the oocyst distribution)
tba_plus <- function(tra, model_fun, mu = mALL, k = rrALL) {
  mf    <- model_fun(j_pos)
  denom <- sum(dnbinom(j_pos, size = k, mu = mu) * mf)
  sapply(tra, function(t) 1 - sum(dnbinom(j_pos, size = k, mu = mu * (1 - t)) * mf) / denom)
}

# NB approximation used inside malariasimulation
calculate_TBA <- function(tra, mx, k) {
  offset <- (k / (k + mx)) ^ k;
  scale <- 1. / (1. - offset);
  tra_transformation <- (k / (k + mx * (1 - tra))) ^ k;
  scale * (tra_transformation - offset)
}

# (m, k) passed to malariasimulation -- keep in sync with the scenario script
# the second ones in the list below are from the dnbinom exact solution using tba_plus() isntead of the rnbinom draws - wrong i think (but extremely similar output to below)
# first ones in list now are those fitted to sims with 20 reps in tbaplus_fitting.R script
model_pars <- tibble(model = model_ids,
                     m = c(mALL, 0.1829639,8.392452e-05),#3.088107e-05, 8.257316e-05),  #6.902878e-06, 7.209415e-05),
                     k = c(rrALL, 0.009642682,9.259975e-06))#1.648756e-06, 9.411515e-06))#4.633258e-07, 1.077945e-05))

tra_fine <- seq(0, 1, 0.01)
tra_pts  <- seq(0.05, 0.95, 0.05)

tba_fit <- bind_rows(lapply(model_ids, function(id) {
  pr <- model_pars[model_pars$model == id, ]
  tibble(model = id, tra = tra_fine, tba = calculate_TBA(tra_fine, pr$m, pr$k))
})) %>% mutate(model = factor(model, levels = lvl))

# tba_pts <- bind_rows(lapply(model_ids, function(id) {
#   tibble(model = id, tra = tra_pts, tba = tba_plus(tra_pts, model_funs[[id]]))
# })) %>% mutate(model = factor(model, levels = lvl))
tba_pts <- readRDS('stan_fitting/outputs/fitting_results.rds') %>%
  filter(type == 'data')

p_c <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'grey40', linewidth = 0.4) +
  geom_line(data = tba_fit, aes(x = tra, y = tba, colour = model), linewidth = 0.7) +
  geom_point(data = tba_pts, aes(x = tra, y = tbaextra, colour = model),
             shape = 21, fill = 'white', size = 1.3, stroke = 0.5) +
  annotate('text', x = 0.30, y = 0.40, label = 'TBA+ = TRA', size = 2.2, colour = 'grey35', hjust = 0.5) +
  annotate('text', x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 2.1, colour = 'grey30', lineheight = 0.95,
           label = 'points: discrete model\nlines: NegBin fit') +
  scale_col_model(model_ids) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = 'Transmission-reducing activity (TRA)',
       y = 'Transmission-blocking activity (TBA+)')

# shared legend, taken from p_c (b and c use the same colour scale)
leg <- cowplot::get_legend(p_c + theme(legend.position = 'bottom'))

p_b_nl <- p_b + theme(legend.position = 'none')
p_c_nl <- p_c + theme(legend.position = 'none')

fig1 <- wrap_plots(p_a, p_b_nl, p_c_nl,
                   wrap_elements(full = leg, ignore_tag = TRUE),
                   design  = 'ABC\n#DD',
                   heights = c(1, 0.1)) +
  plot_annotation(tag_levels = 'a')
# fig1 <- wrap_plots(p_a, p_b, p_c,
#                    design = 'ABC',
#                    heights = c(1, 1.1)) +#p_a,
#   plot_annotation(tag_levels = 'a') +
#   theme(legend.position = 'bottom')
print(fig1)
save_fig(fig1, 'Fig1_oocyst_model', width = 11, height = 5)#width = 7.2, height = 6.4)

# =============================================================================
# FIGURES 2 & 3: EPIDEMIOLOGICAL IMPACT
# =============================================================================
sum_all <- readRDS(input_rds)
if (!is.null(SCENARIO_MATCH)) {
  sum_all <- sum_all %>% filter(grepl(SCENARIO_MATCH, scenario, fixed = TRUE))
}

# ---- Sanity check on stored dose totals (does NOT change any plotted value) -------
# n_tbv_doses is an annual total copied onto each month; if totaldoses was built by summing
# the monthly rows it is ~12x the true total.
doses_annual <- sum_all %>%
  filter(age_group == '0-200') %>%
  distinct(scenario, model, year, n_tbv_doses_med) %>%
  group_by(scenario, model) %>%
  summarise(sum_annual_doses = sum(n_tbv_doses_med), .groups = 'drop')
dose_check <- sum_all %>%
  filter(age_group == '0-200') %>%
  group_by(scenario, model) %>%
  summarise(stored_totaldoses = first(totaldoses_med), .groups = 'drop') %>%
  left_join(doses_annual, by = c('scenario', 'model')) %>%
  mutate(stored_over_annual_sum = stored_totaldoses / sum_annual_doses)
print(as.data.frame(dose_check))
if (any(dose_check$stored_over_annual_sum > 1.5, na.rm = TRUE)) {
  warning('Stored totaldoses is >1.5x the sum of annual doses: per-dose values (Fig 3) look ',
          'inflated in the denominator. Check how totaldoses is built in the scenario script.')
}

# ---- Data preparation ---------------------------------------------------------
prep_epi <- function(sum_all) {
  d <- sum_all %>%
    filter(age_group %in% AGE_GROUPS) %>%
    mutate(age_lab  = factor(as.character(age_group), levels = AGE_GROUPS, labels = AGE_LABELS),
           tbv_ages = factor(tbv_ages, levels = names(TBV_LABELS)),
           tbv_lab  = factor(tbv_ages, levels = names(TBV_LABELS), labels = unname(TBV_LABELS)),
           model    = factor(model, levels = lvl),
           tbv_year = tbv_timing / 365 + 2000) %>%      # primary series; booster = +1 y
    arrange(tbv_ages, model, age_lab, time)
  stopifnot(nrow(d) > 0)
  if (n_distinct(d$scenario) != n_distinct(d$tbv_ages)) {
    stop('More than one scenario per tbv_ages value: tighten SCENARIO_MATCH. Found: ',
         paste(unique(d$scenario), collapse = ', '))
  }

  # No-TBV cumulative incidence (identical in every model row): take it from m0 rows
  base <- d %>%
    filter(model == 'm0') %>%
    transmute(tbv_ages, tbv_lab, age_lab, tbv_year, time,
              model = factor('none', levels = lvl),
              cumul_clinical_med   = cumul_clinical_baseline_med,
              cumul_clinical_lower = cumul_clinical_baseline_lower,
              cumul_clinical_upper = cumul_clinical_baseline_upper,
              clinical_med = clinical_baseline_med,
              clinical_lower = clinical_baseline_lower,
              clinical_upper = clinical_baseline_upper)
  cum_all <- bind_rows(base,
                       d %>% select(tbv_ages, tbv_lab, age_lab, tbv_year, time, model,
                                    clinical_med, clinical_lower, clinical_upper,
                                    cumul_clinical_med, cumul_clinical_lower, cumul_clinical_upper))

  # Range between the three models (per vaccination strategy) = the potential 'top-up' in impact
  env <- d %>%
    group_by(tbv_ages, age_lab, time) %>%
    summarise(lo = min(cumul_clinical_med), hi = max(cumul_clinical_med), .groups = 'drop')

  # Prevalence / infectivity are population-wide: use one age group to avoid duplicates.
  # Models are identical before the vaccine, so that stretch is drawn once as 'No TBV'.
  # (The no-TBV counterfactual after vaccination is not in sum_all.)
  prev <- d %>% filter(age_group == '0-200')
  pre  <- prev %>% filter(model == 'm0', time < tbv_year) %>% mutate(model = factor('none', levels = lvl))
  post <- prev %>% filter(time >= tbv_year - 1 / 12)

  list(d = d, cum_all = cum_all, env = env, prev_df = bind_rows(pre, post),
       tt = distinct(d, tbv_year))
}

# ---- Plot helpers ---------------------------------------------------------------
x_scale    <- scale_x_continuous(breaks = 2002:2005, limits = c(2001.5, 2005), expand = expansion(mult = c(0, 0.01)))
facet_age  <- facet_wrap(~age_lab, ncol = 1, scales = 'free_y')
zero_line  <- geom_hline(yintercept = 0, colour = 'grey55', linewidth = 0.3)
tbv_lines <- function(tt) {
  list(geom_vline(data = tt, aes(xintercept = tbv_year),     linetype = 'dashed', colour = 'grey45', linewidth = 0.3),
       geom_vline(data = tt, aes(xintercept = tbv_year + 1), linetype = 'dotted', colour = 'grey45', linewidth = 0.3))
}
post_tbv <- function(df) filter(df, time >= tbv_year - 1 / 12)
cri_ribbon <- function(data, lo, hi) {
  if (!SHOW_CRI) return(NULL)
  geom_ribbon(data = data,
              aes(x = time, ymin = {{ lo }}, ymax = {{ hi }}, fill = model,
                  group = interaction(model, tbv_ages)),
              alpha = 0.15, colour = NA, inherit.aes = FALSE)
}

# ---- Figure 2 ---------------------------------------------------------------------
make_fig2 <- function(e, title = NULL) {
  tt        <- e$tt
  d         <- post_tbv(e$d)
  cum_post  <- post_tbv(e$cum_all)
  first_tbv <- names(TBV_LABELS)[1]     # no-TBV is the same for every strategy: draw it once

  p_A0 <- ggplot() +
    tbv_lines(tt) +
    geom_line(data = cum_all %>% filter(model !='none'),
              aes(x = time, y = clinical_med, colour = model, linetype = tbv_ages), linewidth = 0.55) +
    geom_line(data = cum_all %>% filter(model == 'none'),
              aes(x = time, y = clinical_med, colour = model), linewidth = 0.55) +
    facet_age + x_scale + scale_col_model() + scale_lt_tbv + #scale_fill_epi +
    labs(x = 'Year', y = 'Monthly incidence\n(cases per 1,000)')

  p_A <- ggplot() +
    tbv_lines(tt) +
    # geom_ribbon(data = e$env,
    #             aes(x = time, ymin = lo, ymax = hi, group = tbv_ages, fill = 'Range across models'),
    #             alpha = 0.25) +
    # cri_ribbon(cum_post, cumul_clinical_lower, cumul_clinical_upper) +
    geom_line(data = filter(cum_post, model != 'none'),
              aes(x = time, y = cumul_clinical_med, colour = model, linetype = tbv_ages), linewidth = 0.55) +
    geom_line(data = filter(e$cum_all, model == 'none', tbv_ages == first_tbv),
              aes(x = time, y = cumul_clinical_med, colour = model), linewidth = 0.55) +
    facet_age + x_scale + scale_col_model() + scale_lt_tbv + #scale_fill_epi +
    labs(x = 'Year', y = 'Cumulative incidence\n(cases per 1,000)')

  p_B <- ggplot(d, aes(x = time, y = ca_clinical_med, colour = model, linetype = tbv_ages)) +
    zero_line + tbv_lines(tt) +
    # cri_ribbon(d, ca_clinical_lower, ca_clinical_upper) +
    geom_line(linewidth = 0.55) +
    facet_age + x_scale + scale_col_model() + scale_lt_tbv + #scale_fill_epi +
    labs(x = 'Year', y = 'Cases averted per 1,000\n(monthly, vs no TBV)')

  p_C <- ggplot(d, aes(x = time, y = cumul_ca_clinical_med, colour = model, linetype = tbv_ages)) +
    zero_line + tbv_lines(tt) +
    # cri_ribbon(d, cumul_ca_clinical_lower, cumul_ca_clinical_upper) +
    geom_line(linewidth = 0.55) +
    facet_age + x_scale + scale_col_model() + scale_lt_tbv + #scale_fill_epi +
    labs(x = 'Year', y = 'Cumulative cases averted\nper 1,000 (vs no TBV)')

  # Prevalence

  p_D <- ggplot(mapping = aes(x = time, colour = model)) +
    tbv_lines(tt) +
    # cri_ribbon(prev_post, lm_prevalence_2_10_lower, lm_prevalence_2_10_upper) +
    # geom_line(data = prev_pre,  aes(y = lm_prevalence_2_10_med), linewidth = 0.55) +
    geom_line(data = e$prev_df, aes(y = lm_prevalence_2_10_med, linetype = tbv_ages), linewidth = 0.55) +
    geom_text(data = tt, aes(x = tbv_year + 0.03, y = Inf), label = 'Primary\nseries',
              hjust = 0, vjust = 1.2, size = 2.1, colour = 'grey30', lineheight = 0.9, inherit.aes = FALSE) +
    geom_text(data = tt, aes(x = tbv_year + 1.03, y = Inf), label = 'Booster',
              hjust = 0, vjust = 1.8, size = 2.1, colour = 'grey30', inherit.aes = FALSE) +
    x_scale + scale_col_model() + scale_lt_tbv + #scale_fill_epi +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = 'Year', y = 'LM PfPR2\u201310')

  # p_E <- ggplot(mapping = aes(x = time, colour = model)) +
  #   tbv_lines(tt) +
  #   # cri_ribbon(prev_post, infectivity_lower, infectivity_upper) +
  #   # geom_line(data = prev_pre,  aes(y = infectivity_med), linewidth = 0.55) +
  #   geom_line(data = e$prev_df, aes(y = infectivity_med, linetype = tbv_ages), linewidth = 0.55) +
  #   x_scale + scale_col_model() + scale_lt_tbv + #scale_fill_epi +
  #   scale_y_continuous(labels = scales::label_number(accuracy = 0.001), limits = c(0, NA),
  #                      expand = expansion(mult = c(0, 0.08))) +
  #   labs(x = 'Year', y = 'Infectivity to mosquitoes')

  design <- 'AABBCCDD\nAABBCCDD\nAABBCCDD\nEEEEEEEE'#'AAAA\nAAAA\nBBBB\nBBBB\nCCCC\nCCCC\nDDDD\nDDDD\nEEEE'#
  wrap_plots(p_A0 + theme(legend.position = 'none'), p_A + theme(legend.position = 'none'),
             p_B + theme(legend.position = 'none'), p_C + theme(legend.position = 'none'),
             p_D + theme(legend.position = 'none'),
             design = design, guides = 'collect') +
    plot_annotation(tag_levels = 'a', title = title) +
    theme(legend.position = 'bottom')
}

# ---- Figure 3 (stored per-dose columns, used as-is) -------------------------------
make_fig3 <- function(e, title = NULL) {
  tt <- e$tt
  d  <- post_tbv(e$d)

  p_a <- ggplot(d, aes(x = time, y = ca_clinical_perdose_med, colour = model, linetype = tbv_ages)) +
    zero_line + tbv_lines(tt) +
    # cri_ribbon(d, ca_clinical_perdose_lower, ca_clinical_perdose_upper) +
    geom_line(linewidth = 0.55) +
    facet_age + x_scale + scale_col_model(model_ids) + scale_lt_tbv + #scale_fill_epi +
    labs(x = 'Year', y = 'Clinical cases averted\nper 1,000 doses (monthly)')

  p_b <- ggplot(d, aes(x = time, y = cumul_ca_clinical_perdose_med, colour = model, linetype = tbv_ages)) +
    zero_line + tbv_lines(tt) +
    # cri_ribbon(d, cumul_ca_clinical_perdose_lower, cumul_ca_clinical_perdose_upper) +
    geom_line(linewidth = 0.55) +
    facet_age + x_scale + scale_col_model(model_ids) + scale_lt_tbv + #scale_fill_epi +
    labs(x = 'Year', y = 'Cumulative clinical cases averted\nper 1,000 doses')

  wrap_plots(p_a, p_b, ncol = 2, guides = 'collect') +
    plot_annotation(tag_levels = 'a', title = title) &
    theme(legend.position = 'bottom')
}

# ---- Summary table ----------------------------------------------------------------
make_table <- function(e) {
  base_pre <- e$cum_all %>% filter(model == 'none', time < tbv_year) %>%
    group_by(tbv_lab, age_lab) %>% summarise(base_pre = max(cumul_clinical_med), .groups = 'drop')
  base_end <- e$cum_all %>% filter(model == 'none') %>%
    group_by(tbv_lab, age_lab) %>% summarise(base_end = max(cumul_clinical_med), .groups = 'drop')

  e$d %>%
    group_by(tbv_lab, model, age_lab) %>%
    slice_max(time, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(base_pre, by = c('tbv_lab', 'age_lab')) %>%
    left_join(base_end, by = c('tbv_lab', 'age_lab')) %>%
    group_by(tbv_lab, age_lab) %>%
    mutate(topup_vs_m0 = cumul_ca_clinical_med - cumul_ca_clinical_med[model == 'm0']) %>%
    ungroup() %>%
    transmute(vaccinated = tbv_lab,
              age_group  = age_lab,
              model      = unname(model_labels[as.character(model)]),
              cases_averted_per_1000        = round(cumul_ca_clinical_med, 1),
              pct_of_post_tbv_cases_averted = round(100 * cumul_ca_clinical_med / (base_end - base_pre), 1),
              topup_vs_m0_per_1000          = round(topup_vs_m0, 1),
              stored_cases_averted_per_1000_doses = round(cumul_ca_clinical_perdose_med, 2))
}

# ---- Run --------------------------------------------------------------------------
e <- prep_epi(sum_all)

fig2 <- make_fig2(e, FIG_TITLE)
print(fig2)
save_fig(fig2, 'Fig2_epi_impact', width = 11, height = 7.6)

fig3 <- make_fig3(e, FIG_TITLE)
print(fig3)
save_fig(fig3, 'Fig3_efficiency', width = 10, height = 5.5)

summary_tbl <- make_table(e)
write.csv(summary_tbl, file.path(fig_dir, 'Table_summary_end_of_horizon.csv'), row.names = FALSE)
print(summary_tbl)
