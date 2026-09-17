
install.packages('pacman')
pacman::p_load(pak, ggplot2, tidyverse, cowplot, devtools, retry, reshape2,
               purrr, ltc)


pak::pak("mrc-ide/site")
pak::pak("mrc-ide/netz")
pak::pak('mrc-ide/malariasimulation@Pfs230_2026')
pak::pak("mrc-ide/postie")
pak::pak("mrc-ide/cali")
