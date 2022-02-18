library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(janitor)
library(lubridate)
library(survival)
library(survminer)

ihq <- read_csv("data_IHQ/IHQ_2022.csv") %>%
  clean_names()

ihq <- ihq %>% mutate(resistant = as.factor(resistant), 
            days = as.numeric(pmin(difftime(last_fllowup, nac_start, units = "day"),
            difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)))

km_fit <- survfit(Surv(days, recurrence) ~ 1, data=ihq)
summary(km_fit, times = c(1,30,60,90*(1:10)))

gsp <- ggsurvplot(km_fit, 
                  data = ihq,
                  palette = '#1997c6',
                  ggtheme = theme_bw(base_size = 16),
                  title = "Relapse-free time",
                  legend = "none") 

png(filename = "plots_IHQ/KM.png", width = 800, height = 400)
print(gsp)
dev.off()

km_resp_fit <- survfit(Surv(days, recurrence) ~ resistant, data=ihq)
gsp <- ggsurvplot(
  km_resp_fit,
  data = ihq,
  palette = c("#b87a7a", "#6d9cbe"),
  conf.int = TRUE,          
  pval = TRUE,             
  risk.table = TRUE, 
  pval.method = TRUE,
  risk.table.col = "strata",
  legend.labs = c("RCB-0/RCB-I", "RCB-II/RCB-III"),
  legend.title = "Response",
  risk.table.height = 0.25, 
  ggtheme = theme_bw(base_size = 14)      
)

png(filename = "plots_IHQ/KM_resistant.png", width = 1000, height = 600)
print(gsp)
dev.off()

ihq <- ihq %>% mutate(gria4_status = if_else(ihc_score_gri_apre_nac >= median(ihq[["ihc_score_gri_apre_nac"]], na.rm = TRUE), "high", "low"),
         ndufaf3_status = if_else(ihc_score_ndu_fpre_nac >=median(ihq[["ihc_score_ndu_fpre_nac"]], na.rm = TRUE), "high", "low"))

km_gria4_fit <- survfit(Surv(days, recurrence) ~ gria4_status, data=ihq)
gsp <- ggsurvplot(
  km_gria4_fit,
  data = ihq,
  palette = c("#b87a7a", "#6d9cbe"),
  conf.int = TRUE,          
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,        
  risk.table.col = "strata",
  legend.labs = c("High", "Low"),
  legend.title = "GRIA4 Status",
  risk.table.height = 0.25, 
  ggtheme = theme_bw(base_size = 14)      
)

png(filename = "plots_IHQ/KM_GRIA4.png", width = 1000, height = 600)
print(gsp)
dev.off()

km_ndufaf3_fit <- survfit(Surv(days, recurrence) ~ ndufaf3_status, data=ihq)
gsp <- ggsurvplot(
  km_ndufaf3_fit,
  data = ihq,
  palette = c("#b87a7a", "#6d9cbe"),
  conf.int = TRUE,          
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,        
  risk.table.col = "strata",
  legend.labs = c("High", "Low"),
  legend.title = "GRIA4 Status",
  risk.table.height = 0.25, 
  ggtheme = theme_bw(base_size = 14)      
)

png(filename = "plots_IHQ/KM_NDUFAF3.png", width = 1000, height = 600)
print(gsp)
dev.off()

