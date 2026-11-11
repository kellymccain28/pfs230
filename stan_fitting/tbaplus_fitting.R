# Script to model the transmission blocking activity considering the transmission reducing activity (TBA+)
# get discrete models then fit m and k parameters to be used in malariasimulation runs

# Bompard parameters for oocyst distribtuion
mALL <- 0.000157   # NB mean
rrALL <- 0.00000495 # NB dispersion (size)
# mm is the mean, rr is k or dispersion paramete (this is eq 8 in 10.1038/s41467-021-21775-3)
x1 <- rnbinom(2*(10^7), mu = mALL, size = rrALL) # get distribution of oocysts with no vaccine (including zeros)

pinf_b0 <- 0.590076 # from b0 here https://github.com/mrc-ide/deterministic-malaria-model/blob/master/R/model_parameters.R


# oocysts --> p(inf) with intermediate steps
# (at the moment, they are not calibrated so that it has the same mean as m0 above)
# (i) Oocysts -> Total SPZ
f2 <- function(x, beta = 10263.54, alpha = 29.95) beta*x / (1 + x/alpha) # default values from stan_fitting/step1_fit_sat.rds (also phi = 1.76)
# (ii) Total SPZ -> SPZ per bite
g3 <- function(x, a = 49.61, b = 4980.93, c = 10143.51) a + b * (x**2/(x**2 + c**2)) # from step1_fit_m2.rds
g6 <- function(x, floor_p = 0, mult = 0.05, scale = 34289.7) x*(floor_p + mult*1/(1 + x/scale) ) # from step2_fit_m1.rds -- not sure if ok
# (iii) SPZ per bite -> Prob(B-S infection)
h <- function(x, baselineprob = 0.01, halfsat = 1000, steep = 1.5) baselineprob + (pinf_b0-baselineprob) * (x**steep/(x**steep + halfsat**steep))
#n=steep;halfsat = K


# step-function model = m0
m0 <- function(x, b0 = pinf_b0) ifelse(x > 0, b0, 0) # bc probability of infection is b0 if at least 1 oocyst
# m1: f2 -> g3 -> h chain
m1 <- function(x) ifelse(x == 0, 0, h(g3(f2(x, beta = 15000), a = 600, b = 1000, c = 20000),
                                      steep = 1.8, baselineprob = 0.05, halfsat = 3000))#, c = 4000

m1_orig <- function(x) ifelse(x==0, 0, h(g3(f2(x))))
# m2: f2 -> g6 -> h chain
m2 <- function(x) ifelse(x == 0, 0, h(g6(f2(x, beta = 35000),
                                         mult = 0.1, scale = 35000),
                                      steep = 7, baselineprob = 0.01, halfsat = 2200))#, mult = 0.3
m2_orig <- function(x) ifelse(x == 0, 0, h(g6(f2(x))))

# get sum of p(inf | o=i) * p(o=i) over all values of the oocyst distribution for each plausible model
x1non0sonly <- x1[x1>0]
m0_sum <- mean(sapply(x1non0sonly, m0))

m1_sum <- mean(sapply(x1non0sonly, m1))
# m1out <- sapply(x1, m1)
# table(m1out)
m2_sum <- mean(sapply(x1non0sonly, m2))

m1_const <- m0_sum / m1_sum # 1.135727
m2_const <- m0_sum / m2_sum # 2.370668

m1_c <- function(x) ifelse(x == 0, 0, h(g3(f2(x, beta = 15000), a = 600, b = 1000, c = 20000),
                                        steep = 1.8, baselineprob = 0.05, halfsat = 3000) * m1_const)
m2_c <- function(x) ifelse(x == 0, 0,  h(g6(f2(x, beta = 35000),
                                            mult = 0.1, scale = 35000),
                                         steep = 7, baselineprob = 0.01, halfsat = 2200) * m2_const)

# # test quickly to find new models
# x <- c(1, 2, 3, 5, 10, 20, 50, 100, 150)
# # plot(data.frame(x, m1 = m1_c(x)))
# # plot(data.frame(x, m2 = m2_c(x)))
# m1_orig = m1_c(x)
# m2_orig = m2_c(x)
#
# #step 1
# totalspz3 <- f2(x, beta = 15000)
# totalspz6 <- f2(x, beta = 35000)
# # step2
# spzperbite3 <- g3(totalspz3, a = 600, b = 1000, c = 20000)
# spzperbite6 <- g6(totalspz6, mult = 0.1, scale = 35000)#20000
#
# # step4
# pinf3 <- h(spzperbite3, steep = 1.8, baselineprob = 0.05, halfsat = 3000)
# pinf6 <- h(spzperbite6, steep = 7, baselineprob = 0.01, halfsat = 2200)
#
# test <- data.frame(oo = x,
#                     totalspz3 = totalspz3,
#                    toalspz6 = totalspz6,
#                     spzperbite3 = spzperbite3,
#                     spzperbite6 =spzperbite6,
#                     pinf3 = pinf3,
#                     pinf6 = pinf6,
#                    b0  = b0,
#                    m1_orig = m1_orig,
#                    m2_orig = m2_orig)
# m1_sum <- mean(test$pinf3)
# m2_sum <- mean(test$pinf6)
#
# m1_const <- m0_sum / m1_sum # 1.135727
# m2_const <- m0_sum / m2_sum # 2.370668
# test$m1 <- test$pinf3 * m1_const
# test$m2 <- test$pinf6 * m2_const
#
# ggplot(test)+
#   geom_line(aes(x = oo, y = totalspz3, color = 'm1')) +
#   geom_line(aes(x = oo, y = totalspz6, color = 'm2')) +
#   # scale_y_log10() +
#   labs(x = 'oocysts', y = 'total spz') +
#   scale_x_continuous(breaks = x)
#
# ggplot(test)+
#   geom_line(aes(x = x, y = spzperbite3, color = 'g3')) +
#   geom_line(aes(x = x, y = spzperbite6, color = 'g6')) +
#   scale_y_log10() +
#   labs(x = 'oocysts', y = 'spz per bite')
#
# ggplot(test)+
#   geom_line(aes(x = x, y = pinf3, color = 'pinf g3')) +
#   geom_line(aes(x = x, y = pinf6, color = 'pinf g6')) +
#   geom_line(aes(x = x, y = b0, color = 'baseline')) +
#   scale_y_log10() +
#   labs(x = 'oocysts', y = 'p inf (non-scaled)') +
#   scale_x_continuous(breaks = x)
#
# ggplot(test)+
#   geom_line(aes(x = x, y = m1, color = 'm1'), linewidth = 1) +
#   geom_line(aes(x = x, y = m2, color = 'm2'), linewidth = 1) +
#   geom_line(aes(x = x, y = m1_orig, color = 'm1'), linetype = 2, linewidth = 1) +
#   geom_line(aes(x = x, y = m2_orig, color = 'm2'), linetype = 2, linewidth = 1) +
#   geom_hline(aes(yintercept = b0, color = 'm0'), linewidth = 1) +
#   scale_y_log10() +
#   labs(x = 'oocysts', y = 'p inf (scaled)', color = 'model') +
#   scale_x_continuous(breaks = x) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme_bw() +
#   theme(legend.position = c(0.7, 0.2))

# Check that p(inf) is never above 1
max(sapply(x1non0sonly, m1_c))
max(sapply(x1non0sonly, m2_c))

# Function to produce 'data' on pr(inf | o = i) for each curve
average_infectivity <- function(tra, model_func, oocysts_counterfactual = x1){

  # rbinom is vectorized over 'size', instead of looping
  o_tbv <- rbinom(length(oocysts_counterfactual),
                  size = oocysts_counterfactual,
                  prob = 1 - tra)

  # #how many data pts
  # ll <- length(oocysts_counterfactual)
  # dy <- data.frame('o' = oocysts_counterfactual)
  # dy$o_tbv <- NA # oocysts that get through the TBV
  # dy$prob_infection <- 0
  #
  # for(i in 1:ll){
  #   dy$o_tbv[i] <- rbinom(1, size = dy$o[i],
  #                         prob = 1-tra)
  #
  # if(dy$o_tbv[i] > 0){
  prob_infection <- model_func(o_tbv)
  # }
  # }
  #print(table(dy$o_tbv == 0))

  return(mean(prob_infection))
}


# Simulate the probability of infection for each number of oocysts with different model functions
store0 <- matrix(0, nrow = 51, ncol = 60)

for(i in 1:51){ # vary TRA
  for(j in 1:20){ # repeat
    store0[i,j] <- average_infectivity(tra = 0.02*(i-1), model_func = m0)
    message(i, ', ', j)
  }

  print(i)
}

gm_m0 <- apply(store0, 1, mean)
str(gm_m0)
gm1_m0 <- apply(store0, 1, quantile, prob = 0.025)
gm2_m0 <- apply(store0, 1, quantile, prob = 0.975)
dg0 <- data.frame('TRA' = 0.02*(seq(1,51,1)-1),
                 'mean' = gm_m0,
                 'ci1' = gm1_m0, 'ci2' = gm2_m0)
dg_m0 <- dg0[dg0$TRA > 0,]
dg_m0$TBA_extra <- 1 - dg_m0$mean / dg0[dg0$TRA==0,]$mean
saveRDS(dg_m0, 'stan_fitting/store_m0_new2.rds')


# x1 <- rnbinom(2*(10^7), mu = mALL, size = rrALL) # get distribution of oocysts with no vaccine (including zeros)

store1 <- matrix(0, nrow = 51, ncol = 30)

for(i in 1:51){ # vary TRA
  for(j in 1:20){ # repeat
    store1[i,j] <- average_infectivity(tra = 0.02*(i-1), model_func = m1_c)
    message(i, ', ', j)
  }

  print(i)
}

gm_m1 <- apply(store1, 1, mean)
str(gm_m1)
gm1_m1 <- apply(store1, 1, quantile, prob = 0.025)
gm2_m1 <- apply(store1, 1, quantile, prob = 0.975)
dg1 <- data.frame('TRA' = 0.02*(seq(1,51,1)-1),
                  'mean' = gm_m1,
                  'ci1' = gm1_m1, 'ci2' = gm2_m1)
dg_m1 <- dg1[dg1$TRA > 0,]
dg_m1$TBA_extra <- 1 - dg_m1$mean / dg1[dg1$TRA==0,]$mean

saveRDS(dg_m1, 'stan_fitting/store_m1_new2.rds')


store2 <- matrix(0, nrow = 51, ncol = 30)
for(i in 1:51){ # vary TRA
  for(j in 1:20){ # repeat
    store2[i,j] <- average_infectivity(tra = 0.02*(i-1), model_func = m2_c)
    message(i, ', ', j)
  }

  print(i)
}

gm_m2 <- apply(store2, 1, mean)
str(gm_m2)
gm1_m2 <- apply(store2, 1, quantile, prob = 0.025)
gm2_m2 <- apply(store2, 1, quantile, prob = 0.975)
dg2 <- data.frame('TRA' = 0.02*(seq(1,51,1)-1),
                  'mean' = gm_m2,
                  'ci1' = gm1_m2, 'ci2' = gm2_m2)
dg_m2 <- dg2[dg2$TRA > 0,]
dg_m2$TBA_extra <- 1 - dg_m2$mean / dg2[dg2$TRA==0,]$mean
saveRDS(dg_m2, 'stan_fitting/store_m2_new2.rds')

# Fitting code
# fitting to the estimated tba plus from the model functions above (this is our "data")
tra_grid <- seq(0.02, 1, 0.02)

calculate_TBA <- function(tra, mx, k) {
  offset <- (k / (k + mx)) ^ k;
  scale <- 1. / (1. - offset);
  tra_transformation <- (k / (k + mx * (1 - tra))) ^ k;
  scale * (tra_transformation - offset)
}
v_effic_ALL <- function(TRA, mx, k
                        ){
  (1/(1-(k/(k+mx))^k)) *
    ( (k/(k+mx*(1-TRA)))^k - (k/(k+mx))^k)
}
# target function
# takes the parameters to fit (mk and k)
# then applies the 'calculate_TBA' function from malariasimulation
dg_m0 <- readRDS('stan_fitting/store_m0_new2.rds')
dg_m1 <- readRDS('stan_fitting/store_m1_new2.rds')
dg_m2 <- readRDS('stan_fitting/store_m2_new2.rds')

tba_plus <- function(tra, model_fun, mu = mALL, k = rrALL) {
  mf    <- model_fun(j_pos)
  denom <- sum(dnbinom(j_pos, size = k, mu = mu) * mf)
  sapply(tra, function(t) 1 - sum(dnbinom(j_pos, size = k, mu = mu * (1 - t)) * mf) / denom)
}


target_fn <- function(par, truth){
  mx_fit <- par[1]
  k_fit  <- par[2]
  pred <- calculate_TBA(tra_grid, mx_fit, k_fit)#calculate_TBA(mx_fit, k_fit, tra_grid)
  sum((pred - truth)^2)
}

truth0 <- dg_m0$TBA_extra#tba_plus(tra_grid, m0)
truth1 <- dg_m1$TBA_extra#tba_plus(tra_grid, m1_c)#
# plot(tra_grid, truth)
fit_mk1 <- optim(par = log(c(46.7, 0.9)), # from https://github.com/mrc-ide/malariasimulation/blob/5f63ccd61c8e893892fd12e7580fc22a7c882285/R/parameters.R#L252
                 fn = function(lp) target_fn(exp(lp), truth = truth1),
                 method = "L-BFGS-B")#,
                 # lower = c(1e-8, 1e-8),
                 # upper = c(Inf, Inf)) # or Nelder-Mead
fit_mk1
mx_fit1 <- exp(fit_mk1$par[1])
k_fit1 <- exp(fit_mk1$par[2])
mx_fit1 #0.3587039
k_fit1 #0.02142727

truth2 <- dg_m2$TBA_extra#tba_plus(tra_grid, m2_c)#
# plot(tra_grid, truth)
fit_mk2 <- optim(par = log(c(mALL, rrALL)),#c(46.7, 0.9),
                 fn = function(lp) target_fn(exp(lp), truth = truth2),
                 method = "L-BFGS-B")#,
                 # lower = c(1e-8, 1e-8),
                 # upper = c(Inf, Inf)) # or Nelder-Mead
fit_mk2
mx_fit2 <- exp(fit_mk2$par[1])
k_fit2 <- exp(fit_mk2$par[2])
mx_fit2 #0.0001570028
k_fit2 #2.049708e-05

pred_m0 <- sapply(tra_grid, calculate_TBA, mALL, rrALL)
pred_m1 <- sapply(tra_grid, calculate_TBA, mx_fit1, k_fit1)
pred_m2 <- sapply(tra_grid, calculate_TBA, mx_fit2, k_fit2)

results <- data.frame(
  tra = tra_grid,
  data_m0 = dg_m0$TBA_extra,#truth0,#
  data_m1 = dg_m1$TBA_extra,#truth1,#
  data_m2 = dg_m2$TBA_extra,#truth2,#
  pred_m0 = pred_m0,
  pred_m1 = pred_m1,
  pred_m2 = pred_m2
) %>%
  pivot_longer(-tra,
               names_to = c('type', 'model'),
               names_pattern = "(.*)_(.*)",
               values_to = 'tbaextra')
saveRDS(results, 'stan_fitting/outputs/fitting_results.rds')

# Plot results (TRA-TBA+) in a nice way :)
tratba <- ggplot(results) +
  geom_line(data = results %>% filter(type == 'pred'),
            aes(x = tra, y = tbaextra, linetype = type, color = model), linewidth = 0.8) +
  geom_abline(slope = 1, intercept= 0, linetype = 'dashed') +
  geom_point(data = results %>% filter(type == 'data'),
             aes(x = tra, y = tbaextra, color = model)) +
  labs(x = 'TRA',
       y = 'TBA+',
       color = "Model",
       linetype = NULL,
       shape = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom')

# Plot m1_c and m2_c curves (prob of BS infection as f(oocyst count))
prinf_oo <- data.frame(
  oocysts = x1non0sonly,
  pinf_m0 = sapply(x1non0sonly, m0),
  pinf_m1 = sapply(x1non0sonly, m1_c),
  pinf_m2 = sapply(x1non0sonly, m2_c)
) %>%
  pivot_longer(-oocysts,
               names_to = 'model',
               values_to = 'pinf')

ooprinf <- ggplot(prinf_oo,
       aes(x = oocysts, y = pinf, color = model, group = model)) +
  geom_line() +
  geom_point() +
  labs(x = 'Oocyst count',
       y = 'P(inf)',
       color = 'Model') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none')

pp <- cowplot::plot_grid(ooprinf, tratba  +
                           theme(legend.position = 'none'))
cowplot::plot_grid(pp, cowplot::get_legend(tratba),
                   ncol = 1,
                   rel_heights = c(1,0.1))
ggsave('stan_fitting/outputs/figure_fittedm1m2_new.png')

