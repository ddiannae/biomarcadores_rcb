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

qpcr <- read_excel("../data/Datos crudos.xlsx", sheet = "qPCR_n=24ER+", 
                   skip = 1) %>% clean_names()

qpcr <- qpcr %>% select(colnames(qpcr)[c(1,35:42)]) 
colnames(qpcr) <- c("folio", "abhd14b", "ndufaf3", "tex264", "ghr", "gria4", 
                    "hgd", "slc12a", "sostcd1")

qpcr <- qpcr %>% inner_join(pacientes, by = "folio") %>%
  mutate(resistance = as.factor(resistance), 
         days = as.numeric(pmin(difftime(last_followup, nac_start, units = "day"),
                                difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)))


## Las variables seleccionadas
qpcr <- qpcr %>% select(recurrence, death, resistance, age, c_t, c_n, tnm, luminal_a, 
                        luminal_b_her2neg, luminal_b_her2pos, folio,
                        ndufaf3, gria4, slc12a, last_followup, nac_start, recurrence_date)

qpcr <- qpcr %>% mutate(resistance = as.factor(resistance), 
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
          gria4_status = as.factor(if_else(gria4 >= median(qpcr[["gria4"]], na.rm = TRUE), "high", "low")),
          slc12a_status = as.factor(if_else(slc12a >=median(qpcr[["slc12a"]], na.rm = TRUE), "high", "low")),
          ndufaf3_status = as.factor(if_else(ndufaf3 >=median(qpcr[["ndufaf3"]], na.rm = TRUE), "high", "low"))) %>%
          select(-ndufaf3, -gria4, -slc12a)


## variables independientes
ivars <- c("recurrence", "death", "resistance")
ovars <- c("days", "last_followup", "nac_start", "recurrence_date", "folio")
## variables dependientes
dvars <- colnames(qpcr)[!colnames(qpcr) %in% c(ivars, ovars)]

recurrence_stats <- map_chr(dvars, ~paste0("Surv(days, recurrence) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = qpcr, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

# A tibble: 10 × 8
# term               estimate std.error robust.se statistic  p.value conf.low conf.high
# <chr>                 <dbl>     <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
# 1 age                   1.49      0.732     0.707     0.564 0.573      0.373      5.95 
# 2 c_t                   3.48      1.07      1.02      1.23  0.220      0.475     25.5  
# 3 c_n                   0.245     1.10      0.441    -3.18  0.00146    0.103      0.583
# 4 tnm                   2.13      1.07      1.00      0.753 0.451      0.299     15.1  
# 5 luminal_a1            7.03      0.744     0.575     3.39  0.000692   2.28      21.7  
# 6 luminal_b_her2neg1    0.360     0.709     0.664    -1.54  0.124      0.0979     1.32 
# 7 luminal_b_her2pos1    0.636     1.07      0.986    -0.459 0.646      0.0921     4.39 
# 8 gria4_statuslow       0.633     0.731     0.708    -0.645 0.519      0.158      2.54 
# 9 slc12a_statuslow      0.714     0.731     0.694    -0.485 0.628      0.183      2.78 
# 10 ndufaf3_statuslow    0.528     0.731     0.695    -0.918 0.359      0.135      2.06 


## Fit del modelo con todas las muestras
cox_model <- coxph(Surv(days, recurrence) ~ c_n+luminal_a+gria4_status+slc12a_status+ndufaf3_status, 
                   data = qpcr,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)

# A tibble: 5 × 8
# term              estimate std.error robust.se statistic p.value conf.low conf.high
# <chr>                <dbl>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
# 1 c_n                  0.206     1.48      1.21     -1.30   0.194    0.0191      2.23
# 2 luminal_a1          19.4       1.32      1.43      2.07   0.0381   1.18      320.  
# 3 gria4_statuslow      2.17      0.946     0.684     1.13   0.256    0.569       8.31
# 4 slc12a_statuslow     0.197     1.08      0.982    -1.66   0.0977   0.0287      1.35
# 5 ndufaf3_statuslow    0.791     1.27      1.27     -0.184  0.854    0.0657      9.53

death_stats <- map_chr(dvars, ~paste0("Surv(days, death) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = qpcr, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

# A tibble: 10 × 8
# term               estimate std.error robust.se statistic   p.value conf.low conf.high
# <chr>                 <dbl>     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
# 1 age                1.45e+ 9  20160.       0.574   36.7    3.05e-295 4.70e+ 8   4.47e+9
# 2 c_t                1.52e+ 0      1.16     1.08     0.386  7.00e-  1 1.82e- 1   1.27e+1
# 3 c_n                1.01e- 1      1.23     0.690   -3.33   8.75e-  4 2.61e- 2   3.89e-1
# 4 tnm                9.26e- 1      1.16     1.07    -0.0717 9.43e-  1 1.14e- 1   7.55e+0
# 5 luminal_a1         1.32e- 8  15630.       0.906  -20.0    3.52e- 89 2.24e- 9   7.80e-8
# 6 luminal_b_her2neg1 1.04e+ 0      1.16     1.10     0.0325 9.74e-  1 1.19e- 1   8.99e+0
# 7 luminal_b_her2pos1 1.46e+ 0      1.16     1.06     0.358  7.20e-  1 1.84e- 1   1.16e+1
# 8 gria4_statuslow    3.81e- 1      1.15     1.14    -0.851  3.95e-  1 4.11e- 2   3.52e+0
# 9 slc12a_statuslow   7.61e-10  20306.       0.584  -35.9    1.29e-282 2.42e-10   2.39e-9
# 10 ndufaf3_statuslow  8.76e- 1      1.00     0.943   -0.140  8.88e-  1 1.38e- 1   5.56e+0

## Hay valores que tienen gran error std.
## Como no salieron variables significativas, ponemos las mismas que en el modelo 
## multivariado anterior.

cox_model_death <- coxph(Surv(days, death) ~ c_n+luminal_a+gria4_status+slc12a_status+ndufaf3_status, 
                   data = qpcr,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)
# # A tibble: 5 × 8
# term                   estimate std.error robust.se statistic   p.value conf.low conf.high
# <chr>                     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
# 1 c_n               0.290              1.23     0.674    -1.84  6.63e-  2 7.73e- 2   1.09e+0
# 2 luminal_a1        0.0000000317   18949.       1.40    -12.3   8.10e- 35 2.03e- 9   4.95e-7
# 3 gria4_statuslow   1.66               1.16     1.13      0.449 6.53e-  1 1.80e- 1   1.53e+1
# 4 slc12a_statuslow  0.00000000300   9877.       0.602   -32.6   6.98e-233 9.22e-10   9.77e-9
# 5 ndufaf3_statuslow 0.548              1.02     0.823    -0.731 4.65e-  1 1.09e- 1   2.75e+0
