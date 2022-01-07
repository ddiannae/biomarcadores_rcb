library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

dir.create("plots_GSE22226")

GPL1708_features <- read_tsv("data_GSE22226/GPL1708_features.tsv")
GPL4133_features <- read_tsv("data_GSE22226/GPL4133_features.tsv")

local_avan <- read_tsv("data_GSE22226/localmente_avanzadas.tsv")

GPL1708_expression <- read_tsv("data_GSE22226/GPL1708_expression.tsv")
GPL4133_expression <- read_tsv("data_GSE22226/GPL4133_expression.tsv")

GPL1708_feature_ids <- GPL1708_features %>% 
  filter(str_detect(gene_symbol, "GRIA4|SLC12A1|NDUFAF3")) %>% 
  select(id, gene_symbol)

GPL4133_feature_ids <- GPL4133_features %>% 
  filter(str_detect(gene_symbol, "GRIA4|SLC12A1|NDUFAF3")) %>% 
  select(id, gene_symbol)

GPL1708_genes <- GPL1708_feature_ids %>% 
  left_join(GPL1708_expression, by = c("id" = "id_ref")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "name", values_to = "expr") %>%
  group_by(gene_symbol, name) %>% 
  summarise(expr = mean(expr))

GPL4133_genes <- GPL4133_feature_ids %>% 
  left_join(GPL4133_expression, by = c("id" = "id_ref")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "name", values_to = "expr") %>%
  group_by(gene_symbol, name) %>% 
  summarise(expr = mean(expr))

expr_genes <- bind_rows(GPL1708_genes, GPL4133_genes) %>% 
  left_join(local_avan %>% select(geo_accession, rs_group), 
            by = c("name" = "geo_accession"))

gp <- ggplot(expr_genes, aes(x = rs_group, y = expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  ylab("Expression") +
  facet_wrap(~gene_symbol, scales = "free_y") + 
  scale_color_brewer(palette = "Set2", name = "", labels = c("Resistentes", "Sensibles"))

png(filename = paste0("plots_GSE22226/genes_expression.png"), width = 800, height = 400)
print(gp)
dev.off()

expr_genes %>% write_tsv("data_GSE22226/GRIA4_SLC12A1.tsv")