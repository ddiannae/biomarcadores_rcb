library(survival)
library(survminer)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(janitor)

qpcr <- read_csv("data_qPCR/QPCR.csv") %>%
  clean_names() 

qpcr <- qpcr %>% mutate(resistant = as.factor(resistant), 
                        days = as.numeric(pmin(difftime(last_fllowup, nac_start, units = "day"),
                                               difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)))

qpcr <- qpcr %>%
  mutate(gria4_status = if_else(gria4 >= median(qpcr[["gria4"]], na.rm = TRUE), "high", "low"),
         slc12a_status = if_else(slc12a >=median(qpcr[["slc12a"]], na.rm = TRUE), "high", "low"),
         ndufaf3_status = if_else(ndufaf3 >=median(qpcr[["ndufaf3"]], na.rm = TRUE), "high", "low"))

qpcr <- qpcr %>% select(folio, days, recurrence, age, rcb_class, c_t, c_n, tnm,
                        her2, ki67, gria4_status, slc12a_status, ndufaf3_status) %>%
  mutate(across(where(is.character), as.factor), 
         across(where(is.factor), as.numeric))

## Fit del modelo con todas las muestras
cox_model <- coxph(Surv(days, recurrence) ~  c_t+gria4_status+slc12a_status+ndufaf3_status+age, 
                   data=qpcr)
summary(cox_model)
# Call:
#   coxph(formula = Surv(days, recurrence) ~ c_t + gria4_status + 
#           slc12a_status + ndufaf3_status + age, data = qpcr)
# 
# n= 38, number of events= 10 
# (4 observations deleted due to missingness)
# 
# coef exp(coef) se(coef)      z Pr(>|z|)   
# c_t             1.32232   3.75210  0.46731  2.830  0.00466 **
#   gria4_status   -0.27833   0.75705  0.87038 -0.320  0.74914   
# slc12a_status   1.95491   7.06328  0.93335  2.095  0.03621 * 
#   ndufaf3_status  0.52964   1.69832  0.92666  0.572  0.56762   
# age             0.15933   1.17272  0.06475  2.461  0.01387 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# c_t                3.752     0.2665    1.5014     9.377
# gria4_status       0.757     1.3209    0.1375     4.169
# slc12a_status      7.063     0.1416    1.1338    44.002
# ndufaf3_status     1.698     0.5888    0.2762    10.442
# age                1.173     0.8527    1.0329     1.331
# 
# Concordance= 0.819  (se = 0.066 )
# Likelihood ratio test= 13.24  on 5 df,   p=0.02
# Wald test            = 9.44  on 5 df,   p=0.09
# Score (logrank) test = 11.52  on 5 df,   p=0.04