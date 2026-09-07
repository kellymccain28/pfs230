# Exploration of Photini's data
library(ggplot2)
library(tidyverse)
library(janitor)

d <- readxl::read_xlsx('data/Pf_Ruptured oocyst_SG spz_final.xlsx') %>%
  clean_names()

names(d)
head(d)

table(d$parasite) # all NF54
range(d$date) # 12 Jan to 3 May 2022
table(d$post_feeding) # day 14 to 20
table(d$feeding)
table(d$exp_id) # 12 experiments
table(d$sample_id, d$exp_id)
table(d$unruptured)
summary(d$sg_spz_q_pcr) # 3 to 2.5e5
table(d$oocyst_count) # unsure what the H and ND and 'bit hard' mean
summary(d$sg_spz_rup_oocyst) # 11 to 15579

dd <- d %>%
  mutate(oocyst_count = as.numeric(oocyst_count))

# Sporozoites per ruptured oocyst
ggplot(data = dd) +
  geom_histogram(aes(x = sg_spz_rup_oocyst,
                     y = after_stat(count / sum(count))))

ggplot(data = dd %>% filter(oocyst_count <=20)) +
  geom_point(aes(x = oocyst_count, y = sg_spz_rup_oocyst)) +
  scale_y_log10()

