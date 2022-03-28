library(survival)
library(survminer)
library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(janitor)
library(purrr)
library(broom)

pacientes <- read_excel("../data/Base de datos-Incanet-Patología-INCAN .xlsx", 
                        skip = 2) %>% clean_names()

ihq <- read_excel("../data/Datos crudos.xlsx", sheet = "IHQ_n=31ER+", skip=3, 
                  col_names = c("folio", "resistant", "NDUFAF3_porcentaje_positividad",
                                "NDUFAF3_positividad", "NDUFAF3_intensidad", "NDUFAF3_score", 
                                "GRIA4_porcentaje_positividad",  "GRIA4_positividad", 
                                "GRIA4_intensidad", "GRIA4_score")) 


ihq <- ihq %>% inner_join(pacientes, by = "folio") %>%
  mutate(resistance = as.factor(resistance), 
         days = as.numeric(pmin(difftime(last_followup, nac_start, units = "day"),
                                difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)))


## Las variables seleccionadas
ihq <- ihq %>% select(recurrence, death, resistance, age, c_t, c_n, tnm, luminal_a, 
                        luminal_b_her2neg, luminal_b_her2pos, folio,
                      NDUFAF3_score, GRIA4_score, last_followup, nac_start, recurrence_date)

ihq <- ihq %>% mutate(resistance = as.factor(resistance), 
                        recurrence = as.factor(recurrence),
                        death = as.factor(death),
                        luminal_a = as.factor(luminal_a),
                        luminal_b_her2neg = as.factor(luminal_b_her2neg),
                        luminal_b_her2pos = as.factor(luminal_b_her2pos),
                        days = as.numeric(pmin(difftime(last_followup, nac_start, units = "day"),
                                               difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)),
                        tnm = as.numeric(grepl("^III", tnm)),
                        c_t = if_else(c_t %in% c(1,2), 0, 1),
                        c_n = if_else(c_n == 0, 0, 1),
                        age = as.numeric(if_else(age < 50, 0, 1)),
                        gria4_status = as.factor(if_else(GRIA4_score >= median(ihq[["GRIA4_score"]], na.rm = TRUE), "high", "low")),
                        ndufaf3_status = as.factor(if_else(NDUFAF3_score >=median(ihq[["NDUFAF3_score"]], na.rm = TRUE), "high", "low"))) %>%
  select(-GRIA4_score, -NDUFAF3_score)


## variables independientes
ivars <- c("recurrence", "death", "resistance")
ovars <- c("days", "last_followup", "nac_start", "recurrence_date", "folio")
## variables dependientes
dvars <- colnames(ihq)[!colnames(ihq) %in% c(ivars, ovars)]

recurrence_stats <- map_chr(dvars, ~paste0("Surv(days, recurrence) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = ihq, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

# # A tibble: 9 × 8
# term             estimate std.error robust.se statistic p.value conf.low conf.high
# <chr>               <dbl>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 age                 1.57      0.764     0.746    0.605    0.545    0.364      6.78
# 2 c_t                 1.05      0.838     0.825    0.0598   0.952    0.208      5.30
# 3 c_n                 0.891     1.08      1.09    -0.106    0.916    0.104      7.60
# 4 tnm                 1.56      1.08      1.10     0.405    0.686    0.182     13.4 
# 5 luminal_a1          0.896     0.837     0.822   -0.134    0.894    0.179      4.49
# 6 luminal_b_her2n…    0.518     0.765     0.726   -0.906    0.365    0.125      2.15
# 7 luminal_b_her2p…    3.63      0.842     0.875    1.47     0.141    0.653     20.2 
# 8 gria4_statuslow     1.05      0.817     0.790    0.0627   0.950    0.223      4.94
# 9 ndufaf3_statusl…    1.19      1.10      0.995    0.172    0.864    0.169      8.34

## Fit del modelo con todas las muestras
cox_model <- coxph(Surv(days, recurrence) ~ c_n+luminal_a+gria4_status+ndufaf3_status, 
                   data = ihq,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)

# # A tibble: 4 × 8
# term             estimate std.error robust.se statistic p.value conf.low conf.high
# <chr>               <dbl>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 c_n                 0.828     1.24       1.75  -0.108     0.914   0.0267     25.7 
# 2 luminal_a1          0.959     0.916      1.00  -0.0414    0.967   0.135       6.82
# 3 gria4_statuslow     1.00      0.882      1.05   0.00203   0.998   0.128       7.87
# 4 ndufaf3_statusl…    1.13      1.15       1.17   0.107     0.915   0.115      11.2 

death_stats <- map_chr(dvars, ~paste0("Surv(days, death) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = ihq, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

# # A tibble: 9 × 8
# term           estimate std.error robust.se statistic   p.value conf.low conf.high
# <chr>             <dbl>     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
#   1 age             1.86e+9  23294.       0.612    34.9   2.66e-266 5.60e+ 8   6.17e+9
# 2 c_t             8.15e-1      1.23     1.19     -0.172 8.64e-  1 7.89e- 2   8.42e+0
# 3 c_n             2.90e-1      1.23     1.23     -1.01  3.14e-  1 2.62e- 2   3.22e+0
# 4 tnm             5.12e-1      1.23     1.24     -0.540 5.89e-  1 4.53e- 2   5.80e+0
# 5 luminal_a1      3.10e-9  15764.       0.659   -29.7   2.35e-194 8.52e-10   1.13e-8
# 6 luminal_b_her…  1.28e+0      1.23     1.18      0.208 8.35e-  1 1.26e- 1   1.30e+1
# 7 luminal_b_her…  5.18e+0      1.23     1.22      1.34  1.79e-  1 4.71e- 1   5.70e+1
# 8 gria4_statusl…  2.05e+0      1.22     1.19      0.601 5.48e-  1 1.97e- 1   2.13e+1
# 9 ndufaf3_statu…  2.82e+0      1.22     1.10      0.942 3.46e-  1 3.26e- 1   2.43e+1

cox_model_death <- coxph(Surv(days, death) ~ c_n+luminal_a+gria4_status+ndufaf3_status, 
                         data = ihq,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)
# # A tibble: 4 × 8
# term            estimate std.error robust.se statistic  p.value conf.low conf.high
# <chr>              <dbl>     <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#   1 c_n             1.49e- 1      1.85     3.20     -0.594 5.52e- 1 2.79e- 4   7.94e+1
# 2 luminal_a1      9.85e-10  16883.       2.46     -8.44  3.16e-17 7.98e-12   1.22e-7
# 3 gria4_statuslow 2.04e+ 0      1.36     0.993     0.718 4.73e- 1 2.91e- 1   1.43e+1
# 4 ndufaf3_status… 5.98e- 1      1.97     3.19     -0.161 8.72e- 1 1.15e- 3   3.12e+2
