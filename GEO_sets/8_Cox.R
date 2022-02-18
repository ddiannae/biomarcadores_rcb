library(survival)
library(survminer)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(janitor)

## Leemos los datos
local_avan <- read_tsv("data/localmente_avanzadas.tsv")

## Convertimos años a días y renombramos el status
local_avan <- local_avan %>% 
  mutate(time = drfs_even_time_years * 365) %>%
  rename("status" = drfs_1_event_0_censored) %>%
  filter(!is.na(rs_group))

## Agregamos la expresión de los genes
rma_data <- read_tsv("data/rs_rma_normalized.tsv")
features <- read_tsv("data/features.tsv")
features <- features %>% filter(str_detect(gene_symbol,
                                           pattern = "GRIA4|NDUFAF3|SLC12A1"))
rma_data <- rma_data %>% filter(feature %in% features$id) %>% 
  inner_join(features %>% select(id, gene_symbol), by = c("feature" = "id")) %>%
  select(gene_symbol, everything(), -feature) %>%
  pivot_longer(-1) %>% pivot_wider(names_from = gene_symbol)

local_avan <- local_avan %>% inner_join(rma_data, by = "name") 
local_avan <- local_avan %>% select(name, time, status, rs_group, age_years, pr_status_ihc, clinical_t_stage, 
                      clinical_ajcc_stage, pathologic_response_rcb_class,
                      pam50_class, GRIA4, NDUFAF3, SLC12A1) %>%
              mutate(across(where(is.character), as.factor), 
                     across(where(is.factor), as.numeric))

## Fit del modelo con todas las muestras
cox_model <- coxph(Surv(time, status) ~ pathologic_response_rcb_class+GRIA4+NDUFAF3+SLC12A1+pam50_class+clinical_ajcc_stage, 
                   data=local_avan)
summary(cox_model)
# n= 169, number of events= 31 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)   
# pathologic_response_rcb_class  1.00256   2.72524  0.30998  3.234  0.00122 **
# GRIA4                         -0.11231   0.89376  0.75715 -0.148  0.88208   
# NDUFAF3                        0.01715   1.01730  0.31161  0.055  0.95611   
# SLC12A1                        0.70341   2.02063  0.68833  1.022  0.30683   
# pam50_class                   -0.32281   0.72411  0.19290 -1.673  0.09424 . 
# clinical_ajcc_stage            0.37303   1.45213  0.19535  1.909  0.05620 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# pathologic_response_rcb_class    2.7252     0.3669    1.4844     5.003
# GRIA4                            0.8938     1.1189    0.2026     3.942
# NDUFAF3                          1.0173     0.9830    0.5523     1.874
# SLC12A1                          2.0206     0.4949    0.5243     7.787
# pam50_class                      0.7241     1.3810    0.4961     1.057
# clinical_ajcc_stage              1.4521     0.6886    0.9902     2.130
# 
# Concordance= 0.735  (se = 0.042 )
# Likelihood ratio test= 20.51  on 6 df,   p=0.002
# Wald test            = 18.47  on 6 df,   p=0.005
# Score (logrank) test = 20.06  on 6 df,   p=0.003

