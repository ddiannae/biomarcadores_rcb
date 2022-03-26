library(dplyr)
library(ggplot2)


load("../data/INCAN_Mama_12-copy.RData")

brca.res.df <- as_tibble(brca.res.df)
brca.res.df <- brca.res.df %>% 
  mutate(regulated = if_else(log2FoldChange > 1 & padj <= 0.05, "Upregulated",
                             if_else(log2FoldChange < 1 & padj <= 0.05, "Downregulated", 
                                     "NS")),
         regulated = as.factor(regulated), 
         label = if_else(external_gene_name %in%  c("SLC12A1", "NDUFAF3", "GRIA4"), 
                         external_gene_name, ""))

regulated.colors <- c("aquamarine4", "black", "brown3")
names(regulated.colors) <- levels(brca.res.df$regulated)

options(ggrepel.max.overlaps = Inf)

volcano <- ggplot(brca.res.df, aes(x = log2FoldChange, y = -log10(padj), 
                        color = regulated, label = label)) +
  geom_point(size = 2) + 
  geom_hline(yintercept = -log10(0.05), linetype= "dashed") +
  geom_vline(xintercept = -1, linetype= "dashed") +
  geom_vline(xintercept = 1, linetype= "dashed") +
  xlim(-10,10) + 
  scale_color_manual(values = regulated.colors) +
  ggrepel::geom_label_repel(color = "black", size = 8) +
  theme_bw(base_size = 24) +
  theme(legend.title = element_blank())


png(filename = "plots/volcano.png", width = 800, height = 400)
print(volcano)
dev.off()

