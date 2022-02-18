library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(ggthemes)

features <- read_tsv("data/features.tsv")
local_avan <- read_tsv("data/localmente_avanzadas.tsv")

rma <- read_tsv("data/rs_rma_normalized.tsv")

feature_ids <- features %>% 
  filter(gene_symbol %in% c("GRIA4", "SLC12A1", "NDUFAF3")) %>% 
  select(id, gene_symbol)

rma_genes <- feature_ids %>% 
  left_join(rma, by = c("id" = "feature")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "name", values_to = "rma_expr") %>% 
  left_join(local_avan %>% select(name, rs_group), by = "name") %>%
  mutate(rs_group = factor(rs_group, levels = c("S", "R")))

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
  labs(y = "RMA normalized", x = "") +
  facet_wrap(~gene_symbol, scales = "free_y") + 
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank()) +
  stat_compare_means(method = "wilcox.test", show.legend = F)+
  scale_color_brewer(palette = "Set2", name = "", labels = c( "Sensitive", "Resistant"))

png(filename = "plots/genes_RMA.png", width = 800, height = 400)
print(gp)
dev.off()

### Combinaciones

## Todos
rma_genes %>% wilcox_test(rma_expr ~ rs_group)
# # A tibble: 1 × 7
# .y.      group1 group2    n1    n2 statistic     p
# * <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>
#   1 rma_expr R      S        429    78     17107 0.752
gp <- ggplot(rma_genes, aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "RMA normalized", x = "", title = "GRIA4, NDUFAF3, SLC12A1 expression") +
  stat_compare_means(method = "wilcox.test", show.legend = F)+
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c( "Sensitive", "Resistant"))

png(filename = "plots/genes_RMA_GRIA4-NDUFAF3-SLC12A1.png", width = 800, height = 400)
print(gp)
dev.off()


## GRIA4 con SLC12A1
gp <- ggplot(rma_genes %>% filter(gene_symbol %in% c("GRIA4", "SLC12A1")), 
             aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "RMA normalized", x = "", title = "GRIA4, SLC12A1 expression") +
  stat_compare_means(method = "wilcox.test", show.legend = F)+
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c( "Sensitive", "Resistant"))

png(filename = "plots/genes_RMA_GRIA4-SLC12A1.png", width = 800, height = 400)
print(gp)
dev.off()

## GRIA4 con NDUFAF3
gp <- ggplot(rma_genes %>% filter(gene_symbol %in% c("GRIA4", "NDUFAF3")), 
             aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "RMA normalized", x = "", title = "GRIA4, NDUFAF3 expression") +
  stat_compare_means(method = "wilcox.test", show.legend = F)+
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots/genes_RMA_GRIA4-NDUFAF3.png", width = 800, height = 400)
print(gp)
dev.off()

## NDUFAF3 con SLC12A1
gp <- ggplot(rma_genes %>% filter(gene_symbol %in% c("NDUFAF3", "SLC12A1")), 
             aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "RMA normalized", x = "", title = "NDUFAF3, SLC12A1 expression") +
  stat_compare_means(method = "wilcox.test", show.legend = F)+
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c( "Sensitive", "Resistant"))

png(filename = "plots/genes_RMA_NDUFAF3-SLC12A1.png", width = 800, height = 400)
print(gp)
dev.off()
