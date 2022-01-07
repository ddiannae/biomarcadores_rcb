library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

features <- read_tsv("data/features.tsv")
local_avan <- read_tsv("data/localmente_avanzadas.tsv")

# local_avan <- read_tsv("data/localmente_avanzadas_pCR_RD.tsv")

gcrma <- read_tsv("data/rs_gcrma_normalized.tsv")
rma <- read_tsv("data/rs_rma_normalized.tsv")

# gcrma <- read_tsv("data/rs_gcrma_normalized_pcr_rd.tsv")
# rma <- read_tsv("data/rs_rma_normalized_pcr_rd.tsv")

feature_ids <- features %>% 
  filter(gene_symbol %in% c("GRIA4", "SLC12A1", "NDUFAF3")) %>% 
  select(id, gene_symbol)

gcrma_genes <- feature_ids %>% 
  left_join(gcrma, by = c("id" = "feature")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "name", values_to = "gcrma_expr") %>% 
  left_join(local_avan %>% select(name, rs_group), by = "name")

rma_genes <- feature_ids %>% 
  left_join(rma, by = c("id" = "feature")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "name", values_to = "rma_expr") %>% 
  left_join(local_avan %>% select(name, rs_group), by = "name")


gp <- ggplot(gcrma_genes, aes(x = rs_group, y = gcrma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  ylab("GCRMA normalized") +
  facet_wrap(~gene_symbol, scales = "free_y") + 
  scale_color_brewer(palette = "Set2", name = "", labels = c("Resistentes", "Sensibles"))

png(filename = paste0("plots/genes_GCRMA.png"), width = 800, height = 400)
# png(filename = paste0("plots/genes_GCRMA_pcr_rd.png"), width = 800, height = 400)
print(gp)
dev.off()

gp <- ggplot(rma_genes, aes(x = rs_group, y = rma_expr, color = rs_group)) +
  geom_boxplot() +
  geom_point() +
  ylab("RMA normalized") +
  facet_wrap(~gene_symbol, scales = "free_y") + 
  scale_color_brewer(palette = "Set2", name = "", labels = c("Resistentes", "Sensibles"))

png(filename = "plots/genes_RMA.png", width = 800, height = 400)
# png(filename = "plots/genes_RMA_pcr_rd.png", width = 800, height = 400)
print(gp)
dev.off()
