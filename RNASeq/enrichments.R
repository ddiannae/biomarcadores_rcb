library(readr)
library(ggthemes)
library(clusterProfiler)
library(janitor)
library(dplyr)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)

load("../data/INCAN_Mama_12-copy.RData")

brca.res.df <- as_tibble(brca.res.df)

deg <- brca.res.df %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, log2FoldChange, padj) %>%
  filter(padj < 0.05)

up <- deg %>% filter(log2FoldChange > 0)
down <- deg %>% filter(log2FoldChange < 0)

universe <-  brca.res.df %>% janitor::clean_names() %>%
  filter(gene_biotype =="protein_coding") %>% as_tibble()

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


png(filename = paste0("plots/KEGG_down_enrichments.png"), width = 800, height = 600)
dotplot(ekk_down, showCategory=30) + 
  ggtitle("KEGG enrichments. Downregulated genes") + 
  theme_bw(base_size = 20)
dev.off()


ekk <- enrichKEGG(gene = deg %>% pull(entrezgene_id),
                  organism = 'hsa',
                  keyType = "ncbi-geneid",
                  universe = universe[["entrezgene_id"]],
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "fdr",
                  qvalueCutoff = 0.2)
head(ekk)

png(filename = paste0("plots/KEGG_enrichments.png"), width = 800, height = 600)
dotplot(ekk, showCategory=30) + 
  ggtitle("KEGG enrichments") + 
  theme_bw(base_size = 20)
dev.off()

ego_up <- enrichGO(gene = up %>% pull(entrezgene_id),
                universe      = universe[["entrezgene_id"]],
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)
head(ego_up)

png(filename = paste0("plots/GO_up_enrichments.png"), width = 800, height = 600)
dotplot(ego_up, showCategory=30) + ggtitle("GO enrichments. Upregulated genes") + 
  theme_bw(base_size = 24)
dev.off()

ego_down <- enrichGO(gene = down %>% pull(entrezgene_id),
                   universe      = universe[["entrezgene_id"]],
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
head(ego_down)

png(filename = paste0("plots/GO_down_enrichments.png"), width = 800, height = 600)
dotplot(ego_down, showCategory=30) + ggtitle("GO enrichments. Downregulated genes") + 
  theme_bw(base_size = 24)
dev.off()

ego <- enrichGO(gene = deg %>% pull(entrezgene_id),
               universe      = universe[["entrezgene_id"]],
               OrgDb         = org.Hs.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2,
               readable      = TRUE)
head(ego)

png(filename = paste0("plots/GO_enrichments.png"), width = 800, height = 800)
dotplot(ego, showCategory=30) + ggtitle("GO enrichments") + 
  theme_bw(base_size = 24)
dev.off()
