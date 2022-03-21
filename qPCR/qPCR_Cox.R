library(survival)
library(survminer)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(janitor)
library(purrr)
library(broom)

qpcr <- read_csv("data_qPCR/originales_qPCR.csv") %>%
  clean_names() 

datos_mama <- read_csv("data_qPCR/base_completa.csv") %>%
  clean_names() 

qpcr <- qpcr %>% inner_join(datos_mama, by = "folio") %>%
  filter(q_pcr32 == 1)

## Son las variables que seleccionaste
qpcr <- qpcr %>% select(folio, recurrence, death, resistant, last_fllowup, 
                        nac_start, recurrence_date,
                        age, c_t, c_n, tnm, er,
                        her2, ndufaf3, gria4, slc12a)

qpcr <- qpcr %>% mutate(resistant = as.factor(resistant), 
          recurrence = as.factor(recurrence),
          death = as.factor(death),
          days = as.numeric(pmin(difftime(last_fllowup, nac_start, units = "day"),
            difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)),
          tnm = as.numeric(grepl("^III", tnm)),
          c_t = if_else(c_t %in% c(1,2), 0, 1),
          c_n = if_else(c_n == 0, 0, 1),
          age = as.numeric(if_else(age < 50, 0, 1)),
          gria4_status = as.factor(if_else(gria4 >= median(qpcr[["gria4"]], na.rm = TRUE), "high", "low")),
          slc12a_status = as.factor(if_else(slc12a >=median(qpcr[["slc12a"]], na.rm = TRUE), "high", "low")),
          ndufaf3_status = as.factor(if_else(ndufaf3 >=median(qpcr[["ndufaf3"]], na.rm = TRUE), "high", "low"))) %>%
          select(-ndufaf3, -gria4, -slc12a)




## variables independientes
ivars <- c("recurrence", "death", "resistant")
ovars <- c("days", "last_fllowup", "nac_start", "recurrence_date", "folio")
## variables dependientes
dvars <- colnames(qpcr)[!colnames(qpcr) %in% c(ivars, ovars)]

recurrence_stats <- map_chr(dvars, ~paste0("Surv(days, recurrence) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = qpcr, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

#   term                   estimate std.error robust.se statistic  p.value conf.low     conf.high
#   <chr>                     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>    <dbl>         <dbl>
# 1 age                       1.64      0.708     0.693    0.713  0.476     4.21e-1         6.38 
# 2 c_t                       5.07      1.06      1.02     1.59   0.112     6.85e-1        37.5  
# 3 c_n                       0.215     1.08      0.404   -3.80   0.000144  9.75e-2         0.475
# 4 tnm                       3.10      1.06      1.01     1.12   0.264     4.26e-1        22.6  
# 5 er                287868280.     9721.        0.500   38.9    0         1.08e+8 767119010.   
# 6 her2                      0.823     0.802     0.770   -0.253  0.800     1.82e-1         3.72 
# 7 gria4_statuslow           0.578     0.731     0.713   -0.769  0.442     1.43e-1         2.34 
# 8 slc12a_statuslow          0.945     0.707     0.683   -0.0829 0.934     2.48e-1         3.60 
# 9 ndufaf3_statuslow         1.18      0.671     0.648    0.258  0.796     3.32e-1         4.21 

### El HR de ER es demasiado grande


## Fit del modelo con todas las muestras
cox_model <- coxph(Surv(days, recurrence) ~ c_n+gria4_status+slc12a_status+ndufaf3_status, 
                   data = qpcr,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)

#   term              estimate std.error robust.se statistic p.value conf.low conf.high
#   <chr>                <dbl>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
# 1 c_n                  0.181     1.30      0.647    -2.64  0.00831   0.0511     0.644
# 2 gria4_statuslow      0.890     0.931     0.815    -0.143 0.886     0.180      4.39 
# 3 slc12a_statuslow     0.681     0.965     0.988    -0.389 0.698     0.0982     4.73 
# 4 ndufaf3_statuslow    0.859     0.921     0.947    -0.160 0.873     0.134      5.50 


death_stats <- map_chr(dvars, ~paste0("Surv(days, death) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = qpcr, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

#   term              estimate std.error robust.se statistic   p.value conf.low conf.high
#   <chr>                <dbl>     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
# 1 age                1.37e+9  20224.       0.559   37.6    3.58e-310  4.59e+8   4.11e+9
# 2 c_t                1.95e+0      1.16     1.10     0.607  5.44e-  1  2.25e-1   1.69e+1
# 3 c_n                7.36e-2      1.23     0.692   -3.77   1.62e-  4  1.90e-2   2.85e-1
# 4 tnm                1.19e+0      1.16     1.09     0.157  8.75e-  1  1.39e-1   1.02e+1
# 5 er                 2.91e+8  14479.       0.619   31.5    2.15e-217  8.64e+7   9.79e+8
# 6 her2               9.76e-1      1.16     1.10    -0.0223 9.82e-  1  1.13e-1   8.43e+0
# 7 gria4_statuslow    3.40e-1      1.15     1.14    -0.947  3.44e-  1  3.65e-2   3.17e+0
# 8 slc12a_statuslow   3.27e-1      1.15     1.12    -1.00   3.16e-  1  3.66e-2   2.91e+0
# 9 ndufaf3_statuslow  1.50e+0      1.00     0.969    0.420  6.74e-  1  2.25e-1   1.00e+1

cox_model_death <- coxph(Surv(days, death) ~ c_n+gria4_status+slc12a_status+ndufaf3_status, 
                   data = qpcr,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)
# A tibble: 4 × 8
#   term              estimate std.error robust.se statistic   p.value conf.low conf.high
#   <chr>                <dbl>     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
# 1 c_n               2.99e-10      1.23     0.684   -32.0   2.61e-225 7.83e-11   1.14e-9
# 2 gria4_statuslow   1.17e+ 0      1.16     1.13      0.140 8.88e-  1 1.27e- 1   1.08e+1
# 3 slc12a_statuslow  2.05e- 9      1.23     0.684   -29.2   7.27e-188 5.36e-10   7.84e-9
# 4 ndufaf3_statuslow 1.10e+ 0      1.02     0.797     0.117 9.07e-  1 2.30e- 1   5.23e+0
