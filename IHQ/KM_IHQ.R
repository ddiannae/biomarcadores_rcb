library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(janitor)
library(lubridate)
library(survival)
library(survminer)

pacientes <- read_excel("../data/Base de datos-Incanet-Patología-INCAN .xlsx", 
                        skip = 2)  %>% clean_names()

ihq <- read_excel("../data/Datos crudos.xlsx", sheet = "IHQ_n=31ER+", skip=3, 
                  col_names = c("folio", "resistant", "NDUFAF3_porcentaje_positividad",
                                "NDUFAF3_positividad", "NDUFAF3_intensidad", "NDUFAF3_score", 
                                "GRIA4_porcentaje_positividad",  "GRIA4_positividad", 
                                "GRIA4_intensidad", "GRIA4_score")) 

ihq <- ihq %>% inner_join(pacientes, by = "folio")
ihq <- ihq %>% mutate(resistant = as.factor(resistant), 
            days = as.numeric(pmin(difftime(last_followup, nac_start, units = "day"),
            difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)))

km_fit <- survfit(Surv(days, recurrence) ~ 1, data=ihq)
summary(km_fit, times = c(1,30,60,90*(1:10)))

gsp <- ggsurvplot(km_fit, 
                  data = ihq,
                  palette = '#1997c6',
                  ggtheme = theme_bw(base_size = 25),
                  title = "Relapse-free time",
                  legend = "none") 

png(filename = "plots_IHQ/KM.png", width = 1000, height = 500)
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
  risk.table.height = 0.30, 
  ggtheme = theme_bw(base_size = 24)      
)

png(filename = "plots_IHQ/KM_resistant_sensitive.png", width = 1000, height = 600)
print(gsp)
dev.off()

ihq <- ihq %>% mutate(gria4_status = if_else(GRIA4_score >= median(ihq[["GRIA4_score"]], na.rm = TRUE), "high", "low"),
         ndufaf3_status = if_else(NDUFAF3_score >=median(ihq[["NDUFAF3_score"]], na.rm = TRUE), "high", "low"))

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
  risk.table.height = 0.30, 
  ggtheme = theme_bw(base_size = 24)      
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
  risk.table.height = 0.30, 
  ggtheme = theme_bw(base_size = 24)      
)

png(filename = "plots_IHQ/KM_NDUFAF3.png", width = 1000, height = 600)
print(gsp)
dev.off()

