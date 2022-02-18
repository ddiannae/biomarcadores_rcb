library(readr)
library(ggthemes)
library(clusterProfiler)
library(janitor)
library(dplyr)
library(pathview)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)

load("data/deg_originales.RData")

deg_originales <- read_tsv("data/299genesDE.csv") %>% janitor::clean_names()

universe <-  brca.res.df %>% janitor::clean_names() %>%
  filter(gene_biotype =="protein_coding") %>% as_tibble()

ekk <- enrichKEGG(gene = deg_originales %>% pull(entrezgene),
                  organism = 'hsa',
                  keyType = "ncbi-geneid",
                  universe = universe[["entrezgene_id"]],
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "fdr",
                  qvalueCutoff = 0.2)
head(ekk)
# Ninguna vía

# geneList <- deg_originales %>% pull(log2fold_change, name = entrezgene)
# pathview(gene.data  = geneList,
#              pathway.id = ekk@result %>% filter(p.adjust  < 0.05) %>% pull(ID),
#              species    = "hsa",
#              limit      = list(gene=max(abs(geneList)), cpd=1),
#              kegg.dir = "kegg_plots")


ego <- enrichGO(gene = deg_originales %>% pull(entrezgene),
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

epa = enrichPathway(deg_originales %>% pull(entrezgene),
                    universe      = universe[["entrezgene_id"]],
                    pvalueCutoff=0.05)

png(filename = paste0("plots/Reactome_enrichments.png"), width = 800, height = 800)
dotplot(epa, showCategory=30, orderBy = "p.adjust", decreasing = FALSE) + 
  ggtitle("Reactome pathways") + 
  theme_bw(base_size = 24)
dev.off()
