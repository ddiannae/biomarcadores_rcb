library(survival)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(survminer)
library(ggthemes)
library(tidyr)

## Leemos los datos
local_avan <- read_tsv("data/localmente_avanzadas.tsv")

## Convertimos años a días y renombramos el status
local_avan <- local_avan %>% 
  mutate(time = drfs_even_time_years * 365) %>%
  rename("status" = drfs_1_event_0_censored) %>%
  filter(!is.na(rs_group))

## Fit del modelo con todas las muestras
km_fit <- survfit(Surv(time, status) ~ 1, data=local_avan)
summary(km_fit, times = c(1,30,60,90*(1:10)))

gsp <- ggsurvplot(km_fit, 
           data = local_avan,
           palette = '#1997c6',
           ggtheme = theme_bw(base_size = 16),
           title = "Distant relapse-free survival",
           legend = "none") 

png(filename = "plots/KM.png", width = 800, height = 400)
print(gsp)
dev.off()

## Fil del modelo separando en grupos Resistente o Sensible
km_resp_fit <- survfit(Surv(time, status) ~ rs_group, data=local_avan)
gsp <- ggsurvplot(
  km_resp_fit,
  data = local_avan,
  palette = c("#b87a7a", "#6d9cbe"),
  conf.int = TRUE,          
  pval = TRUE,             
  risk.table = TRUE, 
  pval.method = TRUE,
  risk.table.col = "strata",
  legend.labs = c("RCB-II/RCB-III", "RCB-0/RCB-I"),
  legend.title = "Response",
  risk.table.height = 0.25, 
  ggtheme = theme_bw(base_size = 14)      
)

png(filename = "plots/KM_RCB.png", width = 1000, height = 600)
print(gsp)
dev.off()

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
local_avan <- local_avan %>%
  mutate(GRIA4_status = if_else(GRIA4 >= median(local_avan[["GRIA4"]]), "high", "low"),
         NDUFAF3_status = if_else(NDUFAF3 >= median(local_avan[["NDUFAF3"]]), "high", "low"),
         SLC12A1_status = if_else(SLC12A1 >=median(local_avan[["SLC12A1"]]), "high", "low"),)


## Fil del modelo separando en high_low
km_GRIA4_fit <- survfit(Surv(time, status) ~ GRIA4_status, data=local_avan)
gsp <- ggsurvplot(
  km_GRIA4_fit,
  data = local_avan,
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

png(filename = "plots/KM_GRIA4.png", width = 1000, height = 600)
print(gsp)
dev.off()

km_NDUFAF3_fit <- survfit(Surv(time, status) ~ NDUFAF3_status, data=local_avan)
gsp <- ggsurvplot(
  km_NDUFAF3_fit,
  data = local_avan,
  palette = c("#b87a7a", "#6d9cbe"),
  conf.int = TRUE,          
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,        
  risk.table.col = "strata",
  legend.labs = c("High", "Low"),
  legend.title = "NDUFAF3 Status",
  risk.table.height = 0.25, 
  ggtheme = theme_bw(base_size = 14)      
)

png(filename = "plots/KM_NDUFAF3.png", width = 1000, height = 600)
print(gsp)
dev.off()

km_SLC12A1_fit <- survfit(Surv(time, status) ~ SLC12A1_status, data=local_avan)
gsp <- ggsurvplot(
  km_SLC12A1_fit,
  data = local_avan,
  palette = c("#b87a7a", "#6d9cbe"),
  conf.int = TRUE,          
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,        
  risk.table.col = "strata",
  legend.labs = c("High", "Low"),
  legend.title = "SLC12A1 Status",
  risk.table.height = 0.25, 
  ggtheme = theme_bw(base_size = 14)      
)

png(filename = "plots/KM_SLC12A1.png", width = 1000, height = 600)
print(gsp)
dev.off()
