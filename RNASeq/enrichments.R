library(readr)
library(clusterProfiler)
library(janitor)
library(dplyr)
library(pathview)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)

load("data_lncrna/deg_originales.RData")

deg <- brca.res.df %>% janitor::clean_names() %>%
  filter(padj <= 0.05 & gene_biotype =="protein_coding") %>% as_tibble()

deg_originales <- read_tsv("data_lncrna/299genesDE.csv") %>% janitor::clean_names()

#deg_originales_protein <- deg_originales %>% filter(gene_biotype =="protein_coding")

deg$ensembl_gene_id %in% deg_originales$ensembl_gene_id

deg <- deg_originales_protein

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

geneList <- deg_originales %>% pull(log2fold_change, name = entrezgene)
pathview(gene.data  = geneList,
             pathway.id = ekk@result %>% filter(p.adjust  < 0.05) %>% pull(ID),
             species    = "hsa",
             limit      = list(gene=max(abs(geneList)), cpd=1),
             kegg.dir = "kegg_plots")



universe <-  brca.res.df %>% janitor::clean_names() %>%
  filter(gene_biotype =="protein_coding") %>% as_tibble()

ego <- enrichGO(gene = deg_originales %>% pull(entrezgene),
               universe      = universe[["entrezgene_id"]],
               OrgDb         = org.Hs.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2,
               readable      = TRUE)
head(ego)

ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(ego, foldChange=geneList) + 
  scale_color_gradient2(name='Log2FC', low='midnightblue', high='red4', 
                        limits = c(-2, 2), oob = scales::squish)
p2 <- cnetplot(ego, categorySize="pvalue", foldChange=geneList) + 
  scale_color_gradient2(name='Log2FC', low='midnightblue', high='red4', 
                        limits = c(-2, 2), oob = scales::squish)
p3 <- cnetplot(ego, foldChange = geneList, circular = TRUE, colorEdge = TRUE) + 
  scale_color_gradient2(name='Log2FC', low='midnightblue', high='red4', 
                        limits = c(-2, 2), oob = scales::squish)

png(filename = "enrich_plots/go_plot.png", width = 1000, height = 400)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
dev.off()

dotplot(ego, showCategory=30) + ggtitle("dotplot for ORA")

epa = enrichPathway(deg %>% pull(entrezgene_id), 
                    pvalueCutoff=0.05)

for(i in epa@result %>% filter(p.adjust  < 0.05) %>% pull(Description)) {
  png(filename = paste0("enrich_plots/", i, "_500_500.png"), width = 500, height = 500)
  print(viewPathway(i, 
              readable = TRUE, 
              foldChange = geneList) + 
    scale_color_gradient2(name='Log2FC', low='midnightblue', high='red4', 
                          limits = c(-2, 2), oob = scales::squish))
  dev.off()
}

