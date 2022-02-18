library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(janitor)

qpcr <- read_csv("data_qPCR/originales_qPCR.csv") %>%
  clean_names() 

datos_mama <- read_csv("data_qPCR/base_completa.csv") %>%
  clean_names() 

qpcr <- qpcr %>% inner_join(datos_mama, by = "folio") %>%
  filter(q_pcr32 == 1)

qpcr <- qpcr %>% select(abhd14b, ndufaf3, tex264, ghr, gria4, 
                        hgd, slc12a, sostcd1, resistant) %>%
  mutate(resistant = as.factor(resistant))
  
dir.create("plots_qPCR")

genes1 <- qpcr %>% select(-resistant, -slc12a, -sostcd1, -tex264) %>% colnames()
genes2 <- c("slc12a", "sostcd1", "tex264") 


lapply(genes1, function(gen) {
  
  gp <- ggplot(qpcr, aes(x = resistant, y = get(gen),
                         color = resistant)) +
    geom_violin(size = 1) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F,) +
    geom_point(position = position_jitter(seed = 1, width = 0.2)) +
    ylab("qPCR expression") +
    xlab("") +
    ggtitle(paste0(toupper(gen), " expression")) +
    theme_bw(base_size = 24) +
    stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                       label.x.npc = "center") +
    scale_x_discrete(labels=c("Sensitive","Resistant")) +
    scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
  
  png(filename = paste0("plots_qPCR/", toupper(gen), ".png"), width = 800, height = 400)
  print(gp)
  dev.off()
  
  gp <- ggplot(qpcr, aes(x = resistant, y = get(gen),
                         color = resistant)) +
    geom_boxplot(size = 1) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F,) +
    geom_point(position = position_jitter(seed = 1, width = 0.2)) +
    ylab("qPCR expression") +
    xlab("") +
    ggtitle(paste0(toupper(gen), " expression")) +
    theme_bw(base_size = 24) +
    stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                       label.x.npc = "center") +
    scale_x_discrete(labels=c("Sensitive","Resistant")) +
    scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
  
  png(filename = paste0("plots_qPCR/", toupper(gen), "_boxplot.png"), width = 800, height = 400)
  print(gp)
  dev.off()
})

lapply(genes2, function(gen) {
  
  gp <- ggplot(qpcr, aes(x = resistant, y = get(gen),
                         color = resistant)) +
    geom_violin(size = 1) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F,) +
    geom_point(position = position_jitter(seed = 1, width = 0.2)) +
    ylab("qPCR expression") +
    xlab("") +
    ggtitle(paste0(toupper(gen), " expression")) +
    geom_signif(comparisons = list(c("0", "1")), 
                map_signif_level=TRUE, size = 1, color = "black", 
                vjust = 0.5, textsize = 7) +
    theme_bw(base_size = 24) +
    theme_bw(base_size = 24) +
    stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                       label.x.npc = "center") +
    scale_x_discrete(labels=c("Sensitive","Resistant")) +
    scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
  
  png(filename = paste0("plots_qPCR/", toupper(gen), ".png"), width = 800, height = 400)
  print(gp)
  dev.off()
  
  gp <- ggplot(qpcr, aes(x = resistant, y = get(gen),
                         color = resistant)) +
    geom_boxplot(size = 1) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F,) +
    geom_point(position = position_jitter(seed = 1, width = 0.2)) +
    ylab("qPCR expression") +
    xlab("") +
    ggtitle(paste0(toupper(gen), " expression")) +
    theme_bw(base_size = 24) +
    geom_signif(comparisons = list(c("0", "1")), 
                map_signif_level=TRUE, size = 1, color = "black", 
                vjust = 0.5, textsize = 7) +
    stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                       label.x.npc = "center") +
    scale_x_discrete(labels=c("Sensitive","Resistant")) +
    scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
  
  png(filename = paste0("plots_qPCR/", toupper(gen), "_boxplot.png"), width = 800, height = 400)
  print(gp)
  dev.off()
  
})
