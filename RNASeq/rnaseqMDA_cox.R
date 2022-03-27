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
                        skip = 2) %>% clean_names() %>%
  mutate(folio = as.character(folio))
  
folios <- read_excel("../data/Datos crudos.xlsx", sheet = "RNAseq_n=66ER+", 
                     range = cell_cols("A:B"), col_names = c("folio", "id")) 
rnaseq <- read_csv("cpm_zcores.tsv") %>% 
  filter(ensembl_id %in% c("ENSG00000178057", "ENSG00000152578", "ENSG00000074803")) %>%
  pivot_longer(cols = starts_with("BCaT"), names_to = "id") %>% 
  pivot_wider(id_cols = id, names_from = ensembl_id, values_from = value) %>%
  rename("ndufaf3" = "ENSG00000178057", "gria4" = "ENSG00000152578", 
         "slc12a1" = "ENSG00000074803") %>%
  inner_join(folios, by ="id")

rnaseq <- pacientes %>% 
  select(folio, recurrence, death, resistance, last_followup, 
         nac_start, recurrence_date, age, c_t, c_n, tnm, 
         luminal_a, luminal_b_her2neg, luminal_b_her2pos) %>%
  inner_join(rnaseq, by ="folio")


rnaseq <- rnaseq %>% mutate(resistance = as.factor(resistance), 
          recurrence = as.factor(recurrence),
          death = as.factor(death),
          days = as.numeric(pmin(difftime(last_followup, nac_start, units = "day"),
            difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)),
          tnm = as.numeric(grepl("^III", tnm)),
          c_t = if_else(c_t %in% c(1,2), 0, 1),
          c_n = if_else(c_n == 0, 0, 1),
          age = as.numeric(if_else(age < 50, 0, 1)),
          gria4_status = as.factor(if_else(gria4 >= median(rnaseq[["gria4"]], na.rm = TRUE), "high", "low")),
          ndufaf3_status = as.factor(if_else(ndufaf3 >=median(rnaseq[["ndufaf3"]], na.rm = TRUE), "high", "low")),
          slc12a_status = as.factor(if_else(slc12a1 >=median(rnaseq[["slc12a1"]], na.rm = TRUE), "high", "low"))) %>%
          select(-ndufaf3, -gria4, -slc12a1, -id)


## variables independientes
ivars <- c("recurrence", "death")
ovars <- c("days", "last_followup", "nac_start", "recurrence_date", "folio")
## variables dependientes
dvars <- colnames(rnaseq)[!colnames(rnaseq) %in% c(ivars, ovars)]

recurrence_stats <- map_chr(dvars, ~paste0("Surv(days, recurrence) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = rnaseq, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

# # A tibble: 11 × 8
# term          estimate std.error robust.se statistic   p.value conf.low conf.high
# <chr>            <dbl>     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
# 1 resistance1    3.57e+0     1.04      1.05     1.21   2.26e-  1  4.54e-1    2.82e1
# 2 age            1.11e+0     0.586     0.589    0.179  8.58e-  1  3.51e-1    3.52e0
# 3 c_t            3.38e+0     1.04      1.04     1.18   2.40e-  1  4.43e-1    2.58e1
# 4 c_n            7.77e+7  7655.        0.501   36.3    6.70e-288  2.91e+7    2.07e8
# 5 tnm            8.52e+7  6362.        0.441   41.4    0          3.59e+7    2.02e8
# 6 luminal_a      1.08e+0     0.586     0.571    0.130  8.97e-  1  3.52e-1    3.29e0
# 7 luminal_b_he…  1.18e+0     0.578     0.567    0.291  7.71e-  1  3.88e-1    3.58e0
# 8 luminal_b_he…  5.48e-1     1.04      1.06    -0.565  5.72e-  1  6.80e-2    4.42e0
# 9 gria4_status…  1.40e+0     0.586     0.588    0.572  5.67e-  1  4.42e-1    4.43e0
# 10 ndufaf3_stat…  1.01e+0     0.578     0.559    0.0229 9.82e-  1  3.39e-1    3.03e0
# 11 slc12a_statu…  7.04e-1     0.586     0.579   -0.607  5.44e-  1  2.26e-1    2.19e0

## Fit del modelo con todas las muestras
cox_model <- coxph(Surv(days, recurrence) ~ c_t+luminal_a+gria4_status+slc12a_status+ndufaf3_status, 
                   data = rnaseq,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)

# # A tibble: 5 × 8
# term             estimate std.error robust.se statistic p.value conf.low conf.high
# <chr>               <dbl>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 c_t                 3.41      1.07      0.977     1.25    0.210    0.502     23.1 
# 2 luminal_a           1.19      0.655     0.550     0.320   0.749    0.406      3.51
# 3 gria4_statuslow     1.36      0.644     0.599     0.508   0.611    0.419      4.38
# 4 slc12a_statuslow    0.560     0.643     0.526    -1.10    0.270    0.200      1.57
# 5 ndufaf3_statusl…    0.931     0.640     0.679    -0.106   0.916    0.246      3.52

death_stats <- map_chr(dvars, ~paste0("Surv(days, death) ~", .x)) %>%
  map(.f = ~coxph(formula = as.formula(.x),  data = rnaseq, id = folio)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE, exponentiate = TRUE)) 

# # A tibble: 11 × 8
# term        estimate std.error robust.se statistic     p.value conf.low conf.high
# <chr>          <dbl>     <dbl>     <dbl>     <dbl>       <dbl>    <dbl>     <dbl>
# 1 resistance1  9.84e-1      1.16     1.16    -0.0137   9.89e-  1 1.01e- 1   9.60e+0
# 2 age          2.33e+0      1.15     1.16     0.729    4.66e-  1 2.40e- 1   2.26e+1
# 3 c_t          2.77e+8  14920.       0.560   34.7      1.72e-264 9.25e+ 7   8.29e+8
# 4 c_n          7.81e+7  13075.       0.644   28.2      2.15e-175 2.21e+ 7   2.76e+8
# 5 tnm          8.60e+7  10874.       0.599   30.5      2.39e-204 2.66e+ 7   2.78e+8
# 6 luminal_a    2.42e-9  12597.       0.533  -37.2      5.95e-303 8.52e-10   6.89e-9
# 7 luminal_b_…  6.88e+8  12240.       0.521   39.0      0         2.48e+ 8   1.91e+9
# 8 luminal_b_…  1.19e-8  11315.       0.605  -30.2      5.85e-200 3.65e- 9   3.90e-8
# 9 gria4_stat…  3.34e-1      1.15     1.15    -0.954    3.40e-  1 3.50e- 2   3.18e+0
# 10 ndufaf3_st…  5.73e+8  12195.       0.524   38.5    Inf.  e-324 2.05e+ 8   1.60e+9
# 11 slc12a_sta…  1.68e-9  12190.       0.524  -38.6      0         6.02e-10   4.69e-9

cox_model_death <- coxph(Surv(days, death) ~c_t+luminal_a+gria4_status+slc12a_status+ndufaf3_status, 
                   data = rnaseq,  id = folio) %>% tidy(conf.int = TRUE, exponentiate = TRUE)

# # A tibble: 5 × 8
# term           estimate std.error robust.se statistic   p.value conf.low conf.high
# <chr>             <dbl>     <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
# 1 c_t             7.85e+8  24792.       0.693   29.5    7.08e-192 2.02e+ 8   3.05e+9
# 2 luminal_a       5.10e-9  18493.       0.868  -22.0    3.18e-107 9.31e-10   2.80e-8
# 3 gria4_statusl…  1.04e+0      1.16     0.966    0.0451 9.64e-  1 1.57e- 1   6.93e+0
# 4 slc12a_status…  1.84e-9  18215.       0.656  -30.7    2.36e-206 5.08e-10   6.66e-9
# 5 ndufaf3_statu…  4.72e+8  17485.       0.658   30.4    1.69e-202 1.30e+ 8   1.71e+9
