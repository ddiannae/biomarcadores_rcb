library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ggpubr)

features <- read_tsv("data_la_extra/features.tsv")
local_avan <- read_tsv("data_la_extra/localmente_avanzadas.tsv")

# local_avan <- read_tsv("data/localmente_avanzadas_pCR_RD.tsv")

rma <- read_tsv("data_la_extra/rs_rma_normalized.tsv")

# rma <- read_tsv("data/rs_rma_normalized_pcr_rd.tsv")

feature_ids <- features %>% 
  filter(gene_symbol %in% c("GRIA4", "SLC12A1", "NDUFAF3")) %>% 
  select(id, gene_symbol)

rma_genes <- feature_ids %>% 
  left_join(rma, by = c("id" = "feature")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "name", values_to = "rma_expr") %>% 
  left_join(local_avan %>% select(name, rs_group), by = "name")

normality_results <- rma_genes %>% group_by(rs_group, gene_symbol) %>% 
  summarise(statistic = shapiro.test(rma_expr)$statistic, 
            value = shapiro.test(rma_expr)$p.value)
## Los datos no son normales

res <- rma_genes %>% group_by(gene_symbol) %>%
  wilcox_test(rma_expr ~ rs_group)
## Localmente avanzadas - No significativos
# <chr>       <chr>    <chr>  <chr>  <int> <int>     <dbl>  <dbl>
# 1 GRIA4       rma_expr R      S        143    26      2002 0.535 
# 2 NDUFAF3     rma_expr R      S        143    26      1636 0.332 
# 3 SLC12A1     rma_expr R      S        143    26      2287 0.0625

#Localmente avanzadas con IIA - No significativos
# gene_symbol .y.      group1 group2    n1    n2 statistic     p
# * <chr>       <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>
#   1 GRIA4       rma_expr R      S        199    46     4987  0.345
# 2 NDUFAF3     rma_expr R      S        199    46     4249  0.45 
# 3 SLC12A1     rma_expr R      S        199    46     5222. 0.137
gp <- ggplot(rma_genes, aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  ylab("RMA normalized") +
  facet_wrap(~gene_symbol, scales = "free_y") + 
  stat_compare_means(method = "wilcox.test")+
  scale_color_brewer(palette = "Set2", name = "", labels = c("Resistentes", "Sensibles"))

png(filename = "plots_la_extra/genes_RMA.png", width = 800, height = 400)
# png(filename = "plots/genes_RMA_pcr_rd.png", width = 800, height = 400)
print(gp)
dev.off()

### Combinaciones

## Todos
rma_genes %>% wilcox_test(rma_expr ~ rs_group)
# # A tibble: 1 × 7
# .y.      group1 group2    n1    n2 statistic     p
# * <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>
#   1 rma_expr R      S        429    78     17107 0.752

## GRIA4 con SLC12A1
gp <- ggplot(rma_genes %>% filter(gene_symbol %in% c("GRIA4", "SLC12A1")), 
             aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "RMA normalized", x = "", title = "GRIA4, SLC12A1") +
  stat_compare_means(method = "wilcox.test")+
  theme_bw(base_size = 18) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Resistentes", "Sensibles"))

png(filename = "plots_la_extra/genes_RMA_GRIA4-SLC12A1.png", width = 800, height = 400)
# png(filename = "plots/genes_RMA_pcr_rd.png", width = 800, height = 400)
print(gp)
dev.off()

## GRIA4 con NDUFAF3
gp <- ggplot(rma_genes %>% filter(gene_symbol %in% c("GRIA4", "NDUFAF3")), 
             aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "RMA normalized", x = "", title = "GRIA4, NDUFAF3") +
  stat_compare_means(method = "wilcox.test")+
  theme_bw(base_size = 18) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Resistentes", "Sensibles"))

png(filename = "plots_la_extra/genes_RMA_GRIA4-NDUFAF3.png", width = 800, height = 400)
print(gp)
dev.off()

## NDUFAF3 con SLC12A1
gp <- ggplot(rma_genes %>% filter(gene_symbol %in% c("NDUFAF3", "SLC12A1")), 
             aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "RMA normalized", x = "", title = "NDUFAF3, SLC12A1") +
  stat_compare_means(method = "wilcox.test")+
  theme_bw(base_size = 18) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Resistentes", "Sensibles"))

png(filename = "plots_la_extra/genes_RMA_NDUFAF3-SLC12A1.png", width = 800, height = 400)
print(gp)
dev.off()
