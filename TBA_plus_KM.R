# Workflow to compare different relationships of P(inf) and oocysts

# P(n | m, k) is the negative binomial distribution of oocyst counts with mean m and dispersion k

# m0 (what is currently in malariasimulation, I think) (Eq 8 in 10.1038/s41467-021-21775-3)
# if a mosquito has oocysts at all, it infects a human with a probability m0
# P(inf | x), where x is oocyst load
# Bompard parameters
mALL <- 0.000157   # NB mean
rrALL <- 0.00000495 # NB dispersion (size)
# mm is the mean, rr is k or dispersion paramete (this is eq 8 in 10.1038/s41467-021-21775-3)
# this is our 'm0', mm shifted down with TRA, then outputs resulting change in probability of being oocyst-pos
# this is a specific case where only if oocysts are entirely killed is there an effect (where the prob is 1 if at least 1 oocyst and 0 if none)
v_effic_ALL <- function(mm = mALL, rr = rrALL,
                        TRA){
  (1/(1-(rr/(rr+mm))^rr)) *
    ( (rr/(rr+mm*(1-TRA)))^rr - (rr/(rr+mm))^rr)
}
v_effic_ALL(TRA = 0.9)
rr = rrALL; mm = mALL
P_baseline = 1 - (rr/(rr+mm))^rr
P_TRA      = 1 - (rr/(rr+mm*(1-TRA)))^rr
Efficacy = 1 - P_TRA/P_baseline
# = (P_baseline - P_TRA)/P_baseline
# = [1 - (rr/(rr+mm))^rr - 1 + (rr/(rr+mm*(1-TRA)))^rr] / [1 - (rr/(rr+mm))^rr]
# = [(rr/(rr+mm*(1-TRA)))^rr - (rr/(rr+mm))^rr] / [1 - (rr/(rr+mm))^rr]


## Get average P(inf) under m0
p0 <- 1 - dnbinom(0, mu = mALL, size = rrALL)
b0 <- 0.590076 # from b0 here https://github.com/mrc-ide/deterministic-malaria-model/blob/master/R/model_parameters.R
target_mean <- b0 * p0 # mean INCLUDING 0s of negbin distribution that the other models should have

# m1 example (at the moment, this is not calibrated so that it has the same mean as m0 above)
# (i) Oocysts -> Total SPZ
f2 <- function(x) 4000*x / (1 + x/25)
# (ii) Total SPZ -> SPZ per bite
g6 <- function(x) x*(0.0001 + 0.14*1/(1 + x/40000) )
# (iii) SPZ per bite -> Prob(B-S infection)
h <- function(x) 0.01 + (b0-0.01) * (x**1.5/(x**1.5 + 1000**1.5))

# f2 -> g6 -> h; needs numerical estimation
x1 <- rnbinom(2*(10^5), mu = mALL, size = rrALL) # get distribution of oocysts with no vaccine (including zeros)
# function to output the iprobability of an infection where mean of this is E[P(inf)|TRA], but we want to get efficacy or TBA+
# which is relative reduction vs no vaccine baseline average_infectivity(tra) → raw mean; TBA+ = 1 - average_infectivity(tra)/average_infectivity(0)
average_infectivity <- function(tra, oocysts_counterfactual = x1){

  dy <- data.frame('o' = oocysts_counterfactual)
  dy$o_tbv <- NA # oocysts that get through the TBV
  dy$spz_tot <- 0
  dy$spz_per_bite <- 0
  dy$prob_infection <- 0

  for(i in 1:length(oocysts_counterfactual)){
    dy$o_tbv[i] <- rbinom(1, size = dy$o[i],
                          prob = 1-tra)
    #
    if(dy$o_tbv[i] > 0){
      dy$spz_tot[i] <- f2(dy$o_tbv[i])
      dy$spz_per_bite[i] <- g6(dy$spz_tot[i])
      dy$prob_infection[i] <- h(dy$spz_per_bite[i])
    }
  }
  return(mean(dy$prob_infection))
}
average_infectivity(tra = 0.2)



# From Claude


v_effic_general <- function(mm, rr, TRA, model_fn, xmax = 2*10^6) {
  xs <- 0:xmax
  pmf0 <- dnbinom(xs, mu = mm,          size = rr)  # baseline oocyst dist.
  pmfT <- dnbinom(xs, mu = mm*(1-TRA),  size = rr)  # post-TRA oocyst dist.

  # prob dist of oocysts times baseline probability of infection
  baseline_mean <- sum(pmf0 * model_fn(xs))
  tra_mean      <- sum(pmfT * model_fn(xs))

  1 - tra_mean / baseline_mean   # TBA+
}

# step-function model = m0 -- should reproduce v_effic_ALL exactly
m0 <- function(x, b0 = 1) ifelse(x > 0, b0, 0) # bc probability of infection is 1 if at least 1 oocyst

# m1: your f2 -> g3 -> h chain
m1 <- function(x) ifelse(x == 0, 0, h(g3(f2(x))))
m2 <- function(x) ifelse(x == 0, 0, h(g6(f2(x))))

v_effic_general(mALL, rrALL, TRA = 0.5, m0)  # sanity check vs v_effic_ALL(TRA=0.5)
v_effic_general(mALL, rrALL, TRA = 0.5, m1_model)    #  new TBA+
v_effic_general(mALL, rrALL, TRA = 0.5, m2)    #  new TBA+

store_m0 <- matrix(0, nrow = 51, ncol = 60)
store_m1 <- matrix(0, nrow = 51, ncol = 60)
store_m2 <- matrix(0, nrow = 51, ncol = 60)
for(i in 1:51){ # vary TRA
  for(j in 1:60){ # repeat
    store_m0[i,j] <- v_effic_general(mALL, rrALL, TRA = 0.02*(i-1), m0)
    store_m1[i,j] <- v_effic_general(mALL, rrALL, TRA = 0.02*(i-1), m1)
    store_m2[i,j] <- v_effic_general(mALL, rrALL, TRA = 0.02*(i-1), m2)
  }

  print(i)
}

gm_m0 <- apply(store_m0, 1, mean)
gm_m1 <- apply(store_m1, 1, mean)
gm_m2 <- apply(store_m2, 1, mean)
dg <- data.frame('TRA' = 0.02*(seq(1,51,1)-1),
                 'TBA_extra_m0' = gm_m0,
                 'TBA_extra_m1' = gm_m1,
                 'TBA_extra_m2' = gm_m2)
ggplot(dg) +
  geom_line(aes(x = TRA, y = mean_m0, color = 'm0')) +
  geom_line(aes(x = TRA, y = mean_m1, color = 'm1')) +
  geom_line(aes(x = TRA, y = mean_m2, color = 'm2')) +
  # geom_line(data = df, aes(x = TRA, y = TBA, color = 'a')) +
  # geom_line(data = df, aes(x = TRA, y = TRA),
  #           color = 'grey55', linetype = 'dotdash') +
  # scale_colour_manual(values = c('black','orchid'),
  #                     labels = c('Current TBA',
  #                                '"Holistic" TBA')) +
  theme_classic() + labs(color = '') +
  theme(legend.position = c(0.18,0.8)) +
  xlab('Transmission Reducing Activity') +
  ylab('Transmission Blocking Activity')

# # Then need to fit the continuous curve -- this doesn't work at the moment
# estimate_TBA_plus <- function(tra, model_fn, mm = mALL, rr = rrALL, xmax = 5000) {
#   xs <- 0:xmax
#   pmf <- dnbinom(xs, mu = mm*(1-tra), size = rr)
#   sum(pmf * model_fn(xs))
# }
# tra_grid <- seq(0, 1, 0.02)
# baseline <- estimate_TBA_plus(0, model_fn)
# tba_plus <- 1 - sapply(tra_grid, estimate_TBA_plus, model_fn = model_fn) / baseline
# fit <- nls(tba_plus ~ 1 - ((1-(rr/(rr+mm))^rr) )^-1 *
#              ((rr/(rr+mm*(1-tra_grid)*a))^rr - (rr/(rr+mm))^rr),
#            start = list(a = 1))
