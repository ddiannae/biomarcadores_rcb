library(readr)
library(dplyr)
library(limma)
library(matrixStats)
library(janitor)

exp_rma <- read_tsv("data/rs_rma_normalized.tsv")
# exp_rma <- read_tsv("data/rs_rma_normalized_pcr_rd.tsv")

genes <- exp_rma %>% pull(feature)
exp_rma <- exp_rma %>% select(-feature) %>% as.matrix()
rownames(exp_rma) <- genes

local_avan <- read_tsv("data/localmente_avanzadas.tsv")
# local_avan <- read_tsv("data/localmente_avanzadas_pCR_RD.tsv")
local_avan <- local_avan %>% filter(name %in% colnames(exp_rma))

medians <- rowMedians(exp_rma)
m_threshold <- 4

png(filename = paste0("plots/medians_RMA.png"), width = 800, height = 400)
# png(filename = paste0("plots/medians_RMA_pcr_rd.png"), width = 800, height = 400)
hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Medians", 
                 border = "gray",
                 xlab = "Median intensities")
abline(v = m_threshold, col = "coral4", lwd = 2)
dev.off()

names(medians) <- genes
ok_medians <- medians[medians > m_threshold]
exp_rma <- exp_rma[names(ok_medians), ]

local_avan <- local_avan %>% mutate(rs_group = as.factor(rs_group))
mm <- model.matrix(~rs_group + 0, data = local_avan)
colnames(mm) <- levels(local_avan$rs_group)
rownames(mm) <- local_avan$name

exp_rma <- exp_rma[, rownames(mm)]

fit <- lmFit(exp_rma, mm)
cm <- makeContrasts(S-R, levels = mm)
fit_RS <- eBayes(contrasts.fit(fit, cm))
tTable <- topTable(fit_RS, number = Inf)

features <- read_tsv("data/features.tsv")
tTable$id <- rownames(tTable)
tTable <- as_tibble(tTable)

tTable <- tTable %>% left_join(features %>% select(id, gene_symbol), 
                               by = "id") %>% select(id, gene_symbol, everything())
tTable %>% clean_names() %>% write_tsv("data/deg_sensibles_vs_resistentes.tsv")
# tTable %>% clean_names() %>% write_tsv("data/deg_sensibles_vs_resistentes_pcr_rd.tsv")

genes <- tTable$gene_symbol
names(genes) <- tTable$id
genes <- genes[rownames(fit_RS$coefficients)]
genes <- ifelse(abs(fit_RS$coefficients) >= 0.5, genes, NA)

png(filename = paste0("plots/volcano_SvsR_top0_5.png"), width = 800, height = 400)
# png(filename = paste0("plots/volcano_SvsR_top0_5_pcr_rd.png"), width = 800, height = 400)
volcanoplot(fit_RS, coef = 1L, style = "p-value", highlight = 100, 
            names = genes,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)                             
dev.off()
