library(dplyr)
library(EnhancedVolcano)

load("../data/INCAN_Mama_12-copy.RData")

brca.res.df <- as_tibble(brca.res.df)

png(filename = "plots/volcano.png", width = 1200, height = 600)
EnhancedVolcano(brca.res.df,
                lab = brca.res.df$external_gene_name,
                selectLab = "",
                x = 'log2FoldChange',
                y = 'padj',
                title ="",
                subtitle = "",
                ylab = "-log10(padj)",
                ylim = c(0,8),
                axisLabSize = 30,
                pCutoff = 0.05,
                pointSize = 5,
                FCcutoff = 2,
                titleLabSize = 30,
                labSize = 12,
                drawConnectors = TRUE,
                boxedLabels = TRUE,
                legendPosition = "top", 
                legendLabSize = 32,
                caption ="",
                col=c('black', 'darkgreen', 'red', 'orange'),
                colAlpha = 4/5)

dev.off()
