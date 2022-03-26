library(sva)
library(edgeR)

load("../data/RNAseqMDA_109.RData")

salmon.counts <- salmon.obj$counts
model.design <- model.matrix(~ E_Receptor + Response, 
                             data = sample.info)

to.keep <- filterByExpr(salmon.counts, 
                        design = model.design,
                        min.count = 2)

### Corregimos el batch del receptor de estrogeno
### esto solo lo hacemos para la visualizacion
### pero no para hacer analisis. 

salmon.counts.corrected <-
  ComBat_seq(salmon.counts, batch = sample.info$E_Receptor,
             group = sample.info$Response)

### Sacamos los CPMs
salmon.counts.corrected.cpm <-
  cpm(salmon.counts.corrected, log = TRUE)

sensitive.samps <- which(sample.info$Response == "sensitive")

### Estandarizamos a zscores
salmon.zscores.cpm <- t(scale(t(salmon.counts.corrected.cpm)))

genes <- c("ENSG00000114779", "ENSG00000178057", "ENSG00000164081",
           "ENSG00000112964", "ENSG00000152578", "ENSG00000113924",
           "ENSG00000171243", "ENSG00000074803")
names(genes) <- c("ABHD14B", "NDUFAF3", "TEX264", "GHR", "GRIA4", "HGD", "SOSTDC1", "SLC12A1") 

tmp.zscores.cpm <- salmon.zscores.cpm[genes,]
rownames(tmp.zscores.cpm) <- names(genes)

all_genes <- rownames(salmon.zscores.cpm)
salmon.zscores.cpm <- as_tibble(salmon.zscores.cpm) 
salmon.zscores.cpm$ensembl_id <- all_genes
write_csv(salmon.zscores.cpm %>% dplyr::select(ensembl_id, everything()), "cpm_zcores.tsv")

library(marray)
library(pheatmap)

fc.lim <- 2 #ceiling(max(abs(tmp.mat)))
breakList <- seq(-fc.lim, fc.lim, by = 0.1)
pal <- maPalette(low="blue", high="red",mid="gray")

ph <- pheatmap(tmp.zscores.cpm,
         annotation_col = sample.info[,c("E_Receptor","Response"), 
                                      drop = F],
         color = colorRampPalette(pal)(length(breakList)),
         breaks = breakList, 
         fontsize = 18)

png(filename = paste0("plots_MDA/EReceptor_corection.png"), width = 1200, 
    height = 600)
ph
dev.off()
