library(readr)
library(dplyr)
library(ggplot2)
library(NOISeq)

local_avan <- read_tsv("data/localmente_avanzadas.tsv")
# local_avan <- read_tsv("data/localmente_avanzadas_pCR_RD.tsv")

exp_rma <- read_tsv("data/rs_rma_normalized.tsv")
# exp_rma <- read_tsv("data/rs_rma_normalized_pcr_rd.tsv")
features <- exp_rma %>% pull(feature)
exp_rma <- exp_rma %>% select(-feature) %>% as.matrix()
rownames(exp_rma) <- features
local_avan <- local_avan %>% filter(name %in% colnames(exp_rma))

PCA <- prcomp(t(exp_rma), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- tibble(name = names(PCA$x[,1]), PC1 = PCA$x[,1], PC2 = PCA$x[,2]) %>%
  left_join(local_avan %>% select(name, rs_group), by = "name")
                     

gp <- ggplot(dataGG, aes(x = PC1, y = PC2, color = rs_group)) +
  geom_point() +
  ggtitle("PCA plot RMA normalization") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_color_brewer(palette = "Set2", name = "", 
                     labels = c("Resistentes", "Sensibles"))

png(filename = paste0("plots/PCA_RMA.png"), width = 800, height = 400)
# png(filename = paste0("plots/PCA_RMA_pcr_rd.png"), width = 800, height = 400)
print(gp)
dev.off()

saveARSyNPlot <- function(expr, samples_data, var) {
 
  mydata <- NOISeq::readData(
    data = expr, 
    factors = samples_data %>% select(!!var))
  
  myARSyN <- ARSyNseq(mydata, norm = "n", logtransf = TRUE)
  
  PCA <- prcomp(t(assayData(myARSyN)$exprs), scale = FALSE)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  dataGG <- tibble(name = names(PCA$x[,1]), PC1 = PCA$x[,1], PC2 = PCA$x[,2]) %>%
    left_join(samples_data %>% select(name, !!var), by = "name")
  
  gp <- ggplot(dataGG, aes(x = PC1, y = PC2, color = get(var))) +
    geom_point() +
    ggtitle("PCA plot RMA normalization") +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    labs(color = var) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = sd_ratio) +
   
  
  png(filename = paste0("plots/PCA_", var, ".png"), width = 800, height = 400)
  print(gp)
  dev.off()
  
  
}

exp_rma <- exp_rma[, local_avan %>% pull(name)]
saveARSyNPlot(exp_rma, local_avan, "pam50_class")
saveARSyNPlot(exp_rma, local_avan, "rs_group")
