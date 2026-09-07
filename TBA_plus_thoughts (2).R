# Script written by JDC for exploration of how the inclusion of downstream effects of a TBV's
# transmission-reducing activity (TRA) will impact efficacy
# previously, only took into account TBA (blocking), so reduction in oocysts did not provide a benefit,
# only full blocking provided a benefit

library(ggplot2)
library(cowplot)

#ZT-NB. Compare to Mali data?

mALL <- 0.000157#Negative binomial parameters
rrALL <- 0.00000495 # #Negative binomial parameters

hist(rnbinom(1000000, mu = 0.000157, size = 0.00000495))
#what proportion @ zero?
dnbinom(0, mu = 0.000157, size = 0.00000495)

#plot ZT

(1/(1-dnbinom(0, mu = 0.000157, size = 0.00000495)))*dnbinom(1:100, mu = 0.000157, size = 0.00000495)
df_NB <- data.frame('o' = seq(1,200),
                    'pmf' = (1/(1-dnbinom(0, mu = 0.000157, size = 0.00000495)))*dnbinom(1:200, mu = 0.000157, size = 0.00000495))
df_NB$cmf <- cumsum(df_NB$pmf)

pl0 <- ggplot(df_NB, aes(x = o, y = pmf)) +
  geom_point() + scale_x_log10() + theme_classic() +
  #scale_y_log10() +
  xlab('Oocyst count') + ylab('Probability')
pl0

pl <- ggplot(df_NB, aes(x = o, y = pmf)) +
  geom_point() + scale_x_log10() +
  scale_y_log10() + theme_classic() +
  xlab('Oocyst count') + ylab('PMF')
pl

plc <- ggplot(df_NB, aes(x = o, y = cmf)) + ylim(c(0,1)) +
  geom_point(colour = 'slateblue') + scale_x_log10() + theme_classic() +
  #scale_y_log10() +
  xlab('Oocyst count') + ylab('Cumulative probability')
plc
cowplot::plot_grid(pl0, plc)
ggsave('Bompard_pmf.pdf', height = 5.2, width = 9)

#plot((1/(1-dnbinom(0, mu = 0.000157, size = 0.00000495)))*dnbinom(1:100, mu = 0.000157, size = 0.00000495))

# v_effic_ALL <- (1/(1-(rrALL/(rrALL+mALL))^rrALL)) *
#   ( (rrALL/(rrALL+mALL*(1-TRA_GSK)))^rrALL - (rrALL/(rrALL+mALL))^rrALL)

# mm is the mean, rr is k or dispersion paramete (this is eq 8 in 10.1038/s41467-021-21775-3)
# this is our 'm0', mm shifted down with TRA, then outputs resulting change in probability of being oocyst-pos
v_effic_ALL <- function(mm = mALL, rr = rrALL,
                        TRA){
  (1/(1-(rr/(rr+mm))^rr)) *
    ( (rr/(rr+mm*(1-TRA)))^rr - (rr/(rr+mm))^rr)
}

#test
v_effic_ALL(TRA = 0.7)
v_effic_ALL(TRA = 0)

#plot(v_effic_ALL(TRA = seq(0,1,0.001)))

df <- data.frame('TRA' = seq(0,1,0.001),
        'TBA' = v_effic_ALL(TRA = seq(0,1,0.001)))
ggplot(df) +
  geom_line(aes(x = TRA, y = TBA), color = 'black') +
  theme_classic() + ylab('Transmission-blocking activity') +
  geom_line(aes(x = TRA, y = TRA), color = 'grey55', linetype = 'dotdash') +
  xlab('Transmission-reducing activity')
ggsave('TRA_TBA_regular.pdf', height = 5.5, width = 6.3)

# (i) Oocysts -> Total SPZ
# (ii) Total SPZ -> SPZ per bite
# (iii) SPZ per bite -> Prob(B-S infection)

# (i)

f <- function(x) 4000*x
f2 <- function(x) 4000*x / (1 + x/25)
d1 <- ggplot() + geom_function(fun = f) +
  geom_function(fun = f2, color = 'slateblue') +
  scale_x_log10(limits = c(1,200)) +
  scale_y_log10() + theme_classic() +
  xlab('Oocyst count') +
  ylab('Total number of sporozoites')
d1
ggsave('step1.pdf', height = 6, width = 7.1)

#check efficiency
f3 <- function(x) f2(x)/x
ins1 <- ggplot() + geom_function(fun = f3, color = 'slateblue') +
  #geom_function(fun = f2, color = 'slateblue') +
  scale_x_log10(limits = c(1,200)) +
  #scale_y_log10() +
  theme_classic() +
  xlab('Oocyst count') +
  ylab('Viable spzs / oocyst')

ggdraw() +
  draw_plot(d1) +
  draw_plot(ins1, x = 0.13, y = 0.58, width = .38, height = .38)
ggsave('step1_ins.pdf', height = 6, width = 7.1)

# (ii)

g0 <- function(x) x # just for y=x line
g <- function(x) x*0.01
g2 <- function(x) 50 + (5000-50) * ((log10(x-100)**3)/((log10(x-100)**3) + (4*3)))
g3 <- function(x) 50 + 4950 * (x**2/(x**2 + 10000**2))
g5 <- function(x) x*(0.00000001 + 0.09*(exp(-x/500000)))
g6 <- function(x) x*(0.0001 + 0.14*1/(1 + x/40000) )
d2 <- ggplot() + geom_function(fun = g) +
  # geom_function(fun = g2, color = 'purple') +
  geom_function(fun = g3, color = 'orangered2') +
  geom_function(fun = g5, color = 'green') +
  geom_function(fun = g6, color = 'blue3') +
  geom_function(fun = g0, color = 'magenta', linetype = 'dashed') +
  scale_x_log10(limits = c(500,1000000)) +
  scale_y_log10() + theme_classic() +
  xlab('Total number of sporozoites') +
  ylab('Number of sporozoites per bite')
d2
ggsave('step2.pdf', height = 6, width = 7.1)

#efficiency for an insert?
g8 <- function(x) 100*g6(x) / x
ins <- ggplot() + geom_function(fun = g8, color = 'blue3') +
  scale_x_log10(limits = c(500,1000000)) +
  #scale_y_log10() +
  theme_classic() +
  xlab('Total number of sporozoites') +
  ylab('Efficiency (%)')

ggdraw() +
  draw_plot(d2) +
  draw_plot(ins, x = 0.11, y = 0.655, width = .378, height = .337)
ggsave('step2_ins.pdf', height = 6, width = 7.1)

# (iii)

h <- function(x) 0.01 + (0.4-0.01) * (x**1.5/(x**1.5 + 1000**1.5))
d3 <- ggplot() +
  geom_function(fun = h, color = 'purple') +
  scale_x_log10(limits = c(1,20000)) +
  #scale_y_log10() +
  theme_classic() +
  xlab('Number of sporozoites per bite') +
  ylab('Prob. of blood-stage infection')
d3
ggsave('step3.pdf', height = 6, width = 7.1)

# Add data points from literature ----
# d1 (oocyst x total spz), d2 (total spz x spz per bite), d3 (spz per bite x p(b-s infection))
d1_literature_spz_oocysts <- readxl::read_xlsx("C:/Users/kem22/OneDrive - Imperial College London/Challenger, Joseph's files - PfVIMT_Joe_Kelly/literature-oocysts-spz.xlsx",
                                                  sheet = 3)
d1_photini <- readxl::read_xlsx('M:/Kelly/postdoc_JoeC/pfs230/data/Pf_Ruptured oocyst_SG spz_final.xlsx') %>%
  clean_names() %>%
  mutate(oocyst_count = as.numeric(oocyst_count))
d2_literature_spz_total_perbite <- readxl::read_xlsx("C:/Users/kem22/OneDrive - Imperial College London/Challenger, Joseph's files - PfVIMT_Joe_Kelly/literature-oocysts-spz.xlsx",
                                                     sheet = 4)
d3_literature_spz_infection <- readxl::read_xlsx("C:/Users/kem22/OneDrive - Imperial College London/Challenger, Joseph's files - PfVIMT_Joe_Kelly/literature-oocysts-spz.xlsx",
                                                 sheet = 5)

d1b <- d1 +
  geom_point(data = d1_literature_spz_oocysts,
             aes(x = oocysts, y = spz_total, color = mosquito)) +
  geom_point(data = d1_photini,
             aes(x = oocyst_count, y = sg_spz_rup_oocyst, color = "Photini's data")) +
  labs(color = NULL) +
  theme(legend.position = c(0.8,0.15))
# d1b
d2b <- d2 +
  geom_point(data = d2_literature_spz_total_perbite %>% filter(mosquito == 'An coluzzi'),
             aes(x = spz_total, y = spz_per_bite, color = reference, shape = mosquito)) +
  theme(legend.position = c(0.2,0.82))

d3b <- d3 +
  geom_point(data = d3_literature_spz_infection,
             aes(x = total_spz, y = pr_infection, color = reference)) +
  theme(legend.position = c(0.2,0.8))
ggsave('step1_literature.pdf', d1b, height = 6, width = 7.1)
ggsave('step2_literature.pdf', d2b, height = 6, width = 7.1)
ggsave('step3_literature.pdf', d3b, height = 6, width = 7.1)


df_NB2 <- df_NB
df_NB2$spz_tot <- f2(df_NB$o)
df_NB2$spz_per_bite <- g6(df_NB2$spz_tot)
df_NB2$prob_infection <- h(df_NB2$spz_per_bite)

i1 <- ggplot(df_NB2) + theme_classic() +
  geom_point(aes(x = spz_tot, y = pmf)) +
  scale_x_log10()
i1

i2 <- ggplot(df_NB2) + theme_classic() +
  geom_point(aes(x = spz_per_bite, y = pmf)) +
  scale_x_log10()
i2

i3 <- ggplot(df_NB2) + theme_classic() +
  geom_point(aes(x = prob_infection, y = pmf)) +
  scale_x_log10()
i3

#expectation of prob (B-S infection)?
sum(df_NB2$pmf * df_NB2$prob_infection)

# Or, do it with samples?

x1 <- rnbinom(2*(10**8), mu = 0.000157, size = 0.00000495)
table(x1>0)
x1 <- x1[x1>0]
table(x1)

dx <- data.frame('o' = x1)
dx$spz_tot <- f2(dx$o)
dx$spz_per_bite <- g6(dx$spz_tot)
dx$prob_infection <- h(dx$spz_per_bite)

head(dx)

j1 <- ggplot(dx) + geom_histogram(aes(x = o,
          y = after_stat(count / sum(count)))) +
  theme_classic() + scale_x_log10() +
  xlab('Oocyst count') + ylab('Relative freq.')
j1

ggdraw() +
  draw_plot(d1) +
  draw_plot(j1, x = 0.14, y = 0.56,
            width = .4, height = .4)

j2 <- ggplot(dx) + geom_histogram(aes(x = spz_tot,
      y = after_stat(count / sum(count)))) +
  theme_classic() + scale_x_log10(limits = c(100,1000000)) +
  xlab('Spz tot') + ylab('Relative freq.')
j2

ggdraw() +
  draw_plot(d2) +
  draw_plot(j2, x = 0.59, y = 0.07,
            width = .37, height = .37)

j3 <- ggplot(dx) + geom_histogram(aes(x = spz_per_bite,
                    y = after_stat(count / sum(count)))) +
  theme_classic() + ylab('Relative freq.') +
  xlab('Spz per bite') + scale_x_log10(limits = c(1,20000)) #
j3

j4 <- ggplot(dx) +
  geom_histogram(aes(x =prob_infection,
          y = after_stat(count / sum(count)))) +
  theme_classic() + ylab('Relative freq.') +
  xlab('Pr (B-S infection)') + xlim(c(0,NA))
j4

ggdraw() +
  draw_plot(d3) +
  draw_plot(j3, x = 0.15, y = 0.5,
            width = .4, height = .4)

mean(dx$prob_infection)


# How to proceed:
# (i) select a TRA;
# (ii) Sample oocyst distribution
# (iii) Note the zeros (or use theoretical work for this?)
# for non-zero oocyst counts, run the multi-step
# procedure. Calculate average Prob(BSI), and compare
# with TBV-free case?
# (iv) reduction in spz due to PEV

#Calculate (counterfactual) oocyst
#counts outside the function?

#Note: my functions don't work for o=0.
# So need to separate those!

average_infectivity <- function(tra, oocysts_counterfactual = x1){

  #how many data pts
  ll <- length(oocysts_counterfactual)
  dy <- data.frame('o' = oocysts_counterfactual)
  dy$o_tbv <- NA # oocysts that get through the TBV
  dy$spz_tot <- 0
  dy$spz_per_bite <- 0
  dy$prob_infection <- 0

  for(i in 1:ll){
    dy$o_tbv[i] <- rbinom(1, size = dy$o[i],
                          prob = 1-tra)
    #
    if(dy$o_tbv[i] > 0){
      dy$spz_tot[i] <- f2(dy$o_tbv[i])
      dy$spz_per_bite[i] <- g3(dy$spz_tot[i])
      dy$prob_infection[i] <- h(dy$spz_per_bite[i])
    }
  }
  #print(table(dy$o_tbv == 0))

  return(mean(dy$prob_infection))
}
average_infectivity(tra = 0.2)

store <- matrix(0, nrow = 51, ncol = 60)

for(i in 1:51){ # vary TRA
  for(j in 1:60){ # repeat
    store[i,j] <- average_infectivity(tra = 0.02*(i-1))
  }

  print(i)
}

gm <- apply(store, 1, mean)
str(gm)
gm1 <- apply(store, 1, quantile, prob = 0.025)
gm2 <- apply(store, 1, quantile, prob = 0.975)
dg <- data.frame('TRA' = 0.02*(seq(1,51,1)-1),
                 'mean' = gm,
                 'ci1' = gm1, 'ci2' = gm2)
ggplot(dg) + geom_line(aes(x = TRA, y = mean)) +
  theme_classic() +
  geom_ribbon(aes(x = TRA, ymin = ci1, ymax = ci2),
    fill = 'slateblue', alpha = .3)

#Use TRA = 0 as a baseline
dg2 <- dg[dg$TRA > 0,]
dg2$TBA_extra <- 1 - dg2$mean / dg[dg$TRA==0,]$mean
ggplot(dg2) + geom_line(aes(x = TRA, y = TBA_extra, color = 'b')) +
  geom_line(data = df, aes(x = TRA, y = TBA, color = 'a')) +
  geom_line(data = df, aes(x = TRA, y = TRA),
            color = 'grey55', linetype = 'dotdash') +
  scale_colour_manual(values = c('black','orchid'),
                      labels = c('Current TBA',
                                 '"Holistic" TBA')) +
  theme_classic() + labs(color = '') +
  theme(legend.position = c(0.18,0.8)) +
  xlab('Transmission Reducing Activity') +
  ylab('Transmission Blocking Activity')
ggsave('TRA_TBA_new.pdf', height = 5.5, width = 6.3)
