library(survival)
library(survminer)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(janitor)
library(purrr)
library(broom)


datos_mama <- read_tsv("data_rnaseq_68/base_completa.csv", skip = 2) %>%
  clean_names() 

datos_mama <- datos_mama %>% filter(rna_seq_68 == 1)

## Son las variables que seleccionaste
datos_mama <- datos_mama %>% select(folio, recurrence, death, resistance, last_followup, 
                        nac_start, recurrence_date,
                        age, c_t, c_n, tnm, er,
                        her2, ndufaf3_rn_aseq, gria4_rn_aseq, slc12a1_rn_aseq) 

datos_mama <- datos_mama %>% mutate(resistance = as.factor(resistance), 
          recurrence = as.factor(recurrence),
          death = as.factor(death),
          days = as.numeric(pmin(difftime(last_followup, nac_start, units = "day"),
            difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)),
          tnm = as.numeric(grepl("^III", tnm)),
          c_t = if_else(c_t %in% c(1,2), 0, 1),
          c_n = if_else(c_n == 0, 0, 1),
          age = as.numeric(if_else(age < 50, 0, 1)),
          gria4_status = as.factor(if_else(gria4_rn_aseq >= median(datos_mama[["gria4_rn_aseq"]], na.rm = TRUE), "high", "low")),
          ndufaf3_status = as.factor(if_else(ndufaf3_rn_aseq >=median(datos_mama[["ndufaf3_rn_aseq"]], na.rm = TRUE), "high", "low")),
          slc12a_status = as.factor(if_else(slc12a1_rn_aseq >=median(datos_mama[["slc12a1_rn_aseq"]], na.rm = TRUE), "high", "low"))) %>%
          select(-ndufaf3_rn_aseq, -gria4_rn_aseq, -slc12a1_rn_aseq)


## variables independientes
ivars <- c("recurrence", "death", "resistance")
ovars <- c("days", "last_followup", "nac_start", "recurrence_date", "folio")
## variables dependientes
dvars <- colnames(datos_mama)[!colnames(datos_mama) %in% c(ivars, ovars)]

recurrence_stats <- map_chr(dvars, ~paste0("Surv(days, recurrence) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = datos_mama, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

# # A tibble: 9 × 8
#   term              estimate std.error robust.se statistic p.value conf.low conf.high
#  <chr>                <dbl>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
# 1 age                  2.18      0.601     0.604     1.29    0.196   0.669       7.13
# 2 c_t                  1.85      0.769     0.758     0.808   0.419   0.418       8.16
# 3 c_n                  1.34      1.04      1.04      0.282   0.778   0.175      10.3 
# 4 tnm                  2.65      1.04      1.04      0.935   0.350   0.344      20.4 
# 5 er                  NA         0        NA        NA      NA      NA          NA   
# 6 her2                 0.464     1.04      1.04     -0.735   0.462   0.0600      3.59
# 7 gria4_statuslow      0.750     0.671     0.657    -0.437   0.662   0.207       2.72
# 8 ndufaf3_statuslow    2.40      0.710     0.666     1.32    0.188   0.652       8.86
# 9 slc12a_statuslow     0.497     0.707     0.695    -1.01    0.315   0.127       1.94

## Fit del modelo con todas las muestras
cox_model <- coxph(Surv(days, recurrence) ~ c_t+c_n+gria4_status+slc12a_status+ndufaf3_status, 
                   data = datos_mama,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)

# # A tibble: 5 × 8
#   term                   estimate std.error robust.se statistic   p.value     conf.low    conf.high
#   <chr>                     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>        <dbl>        <dbl>
# 1 c_t                       3.22      1.07      1.03     1.14   2.56e-  1        0.428        24.2 
# 2 c_n               128706132.     8439.        0.611   30.6    5.36e-205 38846372.    426430256.  
# 3 gria4_statuslow           1.06      0.743     0.705    0.0891 9.29e-  1        0.267         4.24
# 4 slc12a_statuslow          0.485     0.746     0.714   -1.01   3.11e-  1        0.120         1.97
# 5 ndufaf3_statuslow         2.28      0.752     0.708    1.16   2.44e-  1        0.569         9.12

death_stats <- map_chr(dvars, ~paste0("Surv(days, death) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = datos_mama, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

# # A tibble: 9 × 8
#   term              estimate std.error robust.se statistic    p.value  conf.low conf.high
#   <chr>                <dbl>     <dbl>     <dbl>     <dbl>      <dbl>     <dbl>     <dbl>
# 1 age                4.78e+0      1.10     1.10      1.42   1.55e-  1  5.54e- 1   4.12e+1
# 2 c_t                1.72e+0      1.10     1.08      0.503  6.15e-  1  2.09e- 1   1.42e+1
# 3 c_n                5.83e-1      1.10     1.06     -0.506  6.13e-  1  7.25e- 2   4.70e+0
# 4 tnm                1.12e+0      1.10     1.08      0.106  9.16e-  1  1.36e- 1   9.27e+0
# 5 er                NA            0       NA        NA     NA         NA         NA      
# 6 her2               1.12e+0      1.10     1.07      0.104  9.17e-  1  1.36e- 1   9.18e+0
# 7 gria4_statuslow    4.67e-1      1.22     1.20     -0.633  5.27e-  1  4.44e- 2   4.92e+0
# 8 ndufaf3_statuslow  6.58e+8  14106.       0.598    34.0    1.13e-252  2.04e+ 8   2.13e+9
# 9 slc12a_statuslow   1.70e-9  14076.       0.599   -33.7    8.88e-249  5.26e-10   5.51e-9

cox_model_death <- coxph(Surv(days, death) ~ c_t+c_n+gria4_status+slc12a_status+ndufaf3_status, 
                   data = datos_mama,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)
# # A tibble: 5 × 8
#   term              estimate std.error robust.se statistic   p.value conf.low conf.high
#   <chr>                <dbl>     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
# 1 c_t               1.48e+ 9  29319.       0.702    30.1   1.24e-198 3.74e+ 8  5.87e+ 9
# 2 c_n               2.65e+ 9  44516.       0.786    27.6   1.53e-167 5.67e+ 8  1.24e+10
# 3 gria4_statuslow   1.48e+ 0      1.23     0.957     0.411 6.81e-  1 2.27e- 1  9.66e+ 0
# 4 slc12a_statuslow  6.89e-10  21495.       0.644   -32.7   3.07e-235 1.95e-10  2.44e- 9
# 5 ndufaf3_statuslow 1.28e+ 9  20732.       0.653    32.1   3.42e-226 3.57e+ 8  4.62e+ 9