library(readr)
library(ggthemes)
library(clusterProfiler)
library(janitor)
library(dplyr)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)
library(purrr)

load("../data/INCAN_Mama_12-copy.RData")

brca.res.df <- as_tibble(brca.res.df)

deg <- brca.res.df %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, log2FoldChange, padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2)

up <- deg %>% filter(log2FoldChange > 0)
down <- deg %>% filter(log2FoldChange < 0)

universe <-  brca.res.df %>% janitor::clean_names() %>%
   as_tibble()

ekk_up <- enrichKEGG(gene = up %>% pull(entrezgene_id),
                  organism = 'hsa',
                  keyType = "ncbi-geneid",
                  universe = universe[["entrezgene_id"]],
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "fdr",
                  qvalueCutoff = 0.2)
head(ekk_up)
## Ninguna via 

ekk_down <- enrichKEGG(gene = down %>% pull(entrezgene_id),
                     organism = 'hsa',
                     keyType = "ncbi-geneid",
                     universe = universe[["entrezgene_id"]],
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "fdr",
                     qvalueCutoff = 0.2)
head(ekk_down)


ego_up <- enrichGO(gene = up %>% pull(entrezgene_id),
                universe      = universe[["entrezgene_id"]],
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)
head(ego_up)

ego_down <- enrichGO(gene = down %>% pull(entrezgene_id),
                   universe      = universe[["entrezgene_id"]],
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
head(ego_down)

go_results <- ego_up@result %>% filter(p.adjust < 0.05) %>%
  mutate(category = "Upregulated") %>% bind_rows(
ego_down@result %>% filter(p.adjust < 0.05) %>%
  mutate(category = "Downregulated")) %>% 
  clean_names() %>%
  mutate(order_label =  paste(category, description, sep = "-"))

go_results$gr <- map_dbl(go_results$gene_ratio, function(x){
  nums <- as.numeric(unlist(strsplit(x, "/")))
  return(nums[1]/nums[2])
})

png(filename = paste0("plots/GO_enrichments.png"), width = 950, height = 600)
ggplot(go_results, aes(x = category, y = order_label, color = p_adjust))  +
  geom_point(aes(size = count)) +
  theme_bw(base_size = 26) +
  theme(legend.key.size  = unit(0.8, 'cm'))+
  scale_y_discrete(labels = ~ map_chr(., ~ unlist(strsplit(., "-"))[2]), 
                   limits=rev) +
  labs(x = "", y = element_blank(), title = "Gene Ontology enrichments") + 
  scale_color_viridis_c(name = "Padj", direction = -1, option = "plasma",
                        breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), 
                        limits = c(0, 0.05), 
                        labels = c("", 0.01, 0.02, 0.03, 0.04,  0.05)) +
  scale_size(name = "Count",range = c(1, 10), breaks = c(1, 3, 5, 7, 9, 11),
             limits = c(1, 11))
dev.off()


kegg_results <- ekk_up@result %>% filter(p.adjust < 0.05) %>%
  mutate(category = "Upregulated") %>% bind_rows(
    ekk_down@result %>% filter(p.adjust < 0.05) %>%
      mutate(category = "Downregulated")) %>% 
  clean_names() %>%
  mutate(order_label =  paste(category, description, sep = "-"))
  
kegg_results$gr <- map_dbl(kegg_results$gene_ratio, function(x){
  nums <- as.numeric(unlist(strsplit(x, "/")))
  return(nums[1]/nums[2])
})

png(filename = paste0("plots/KEGG_enrichments.png"), width = 950, height = 500)
ggplot(kegg_results, aes(x = category, y = order_label, color = p_adjust))  +
  geom_point(aes(size = count)) +
  theme_bw(base_size = 26) +
  theme(legend.key.size  = unit(0.8, 'cm'))+
  scale_y_discrete(labels = ~ map_chr(., ~ unlist(strsplit(., "-"))[2]), 
                   limits=rev) +
  labs(x = "", y = element_blank(), title = "KEGG enrichments") + 
  scale_color_viridis_c(name = "Padj", direction = -1, option = "plasma",
                        breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), 
                        limits = c(0, 0.05), 
                        labels = c("", 0.01, 0.02, 0.03, 0.04,  0.05)) +
  scale_size(name = "Count",range = c(1, 10), breaks = c(1, 3, 5, 7, 9, 11),
             limits = c(1, 11))
dev.off()
