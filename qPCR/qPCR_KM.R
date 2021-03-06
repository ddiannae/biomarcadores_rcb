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

km_fit <- survfit(Surv(days, recurrence) ~ 1, data=qpcr)
summary(km_fit, times = c(1,30,60,90*(1:10)))

gsp <- ggsurvplot(km_fit, 
                  data = qpcr,
                  palette = '#1997c6',
                  ggtheme = theme_bw(base_size = 24),
                  title = "Relapse-free time",
                  legend = "none") 


png(filename = "plots_qPCR/KM.png", width = 800, height = 400)
print(gsp)
dev.off()

## Fil del modelo separando en grupos Resistente o Sensible
km_resp_fit <- survfit(Surv(days, recurrence) ~ resistance, data=qpcr)
gsp <- ggsurvplot(
  km_resp_fit,
  data = qpcr,
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

png(filename = "plots_qPCR/KM_resistant.png", width = 1000, height = 600)
print(gsp)
dev.off()

qpcr <- qpcr %>%
  mutate(gria4_status = if_else(gria4 >= median(qpcr[["gria4"]], na.rm = TRUE), "high", "low"),
         slc12a_status = if_else(slc12a >=median(qpcr[["slc12a"]], na.rm = TRUE), "high", "low"),
         ndufaf3_status = if_else(ndufaf3 >=median(qpcr[["ndufaf3"]], na.rm = TRUE), "high", "low"))

km_gria4_fit <- survfit(Surv(days, recurrence) ~ gria4_status, data=qpcr)
gsp <- ggsurvplot(
  km_gria4_fit,
  data = qpcr,
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

png(filename = "plots_qPCR/KM_GRIA4.png", width = 1000, height = 600)
print(gsp)
dev.off()

km_slc12a_fit <- survfit(Surv(days, recurrence) ~ slc12a_status, data=qpcr)
gsp <- ggsurvplot(
  km_slc12a_fit,
  data = qpcr,
  palette = c("#b87a7a", "#6d9cbe"),
  conf.int = TRUE,          
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,        
  risk.table.col = "strata",
  legend.labs = c("High", "Low"),
  legend.title = "SLC12A Status",
  risk.table.height = 0.30, 
  ggtheme = theme_bw(base_size = 24)      
)

png(filename = "plots_qPCR/KM_SLC12A.png", width = 1000, height = 600)
print(gsp)
dev.off()

km_ndufaf3_fit <- survfit(Surv(days, recurrence) ~ ndufaf3_status, data=qpcr)
gsp <- ggsurvplot(
  km_ndufaf3_fit,
  data = qpcr,
  palette = c("#b87a7a", "#6d9cbe"),
  conf.int = TRUE,          
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,        
  risk.table.col = "strata",
  legend.labs = c("High", "Low"),
  legend.title = "NDUFAF3 Status",  
  risk.table.height = 0.30, 
  ggtheme = theme_bw(base_size = 24)      
)

png(filename = "plots_qPCR/KM_NDUFAF3.png", width = 1000, height = 600)
print(gsp)
dev.off()
