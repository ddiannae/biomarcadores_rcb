library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(janitor)
library(tidyr)

pacientes <- read_excel("../data/Base de datos-Incanet-Patología-INCAN .xlsx", 
                        skip = 2) %>% clean_names()

qpcr <- read_excel("../data/Datos crudos.xlsx", sheet = "qPCR_n=24ER+", 
          skip = 1) %>% clean_names()

qpcr <- qpcr %>% select(colnames(qpcr)[c(1,43:50)]) 
colnames(qpcr) <- c("folio", "abhd14b", "ndufaf3", "tex264", "ghr", "gria4", 
                    "hgd", "slc12a", "sostcd1")

qpcr <- qpcr %>% inner_join(pacientes, by = "folio")
qpcr <- qpcr %>% select(abhd14b, ndufaf3, tex264, ghr, gria4, 
                        hgd, slc12a, sostcd1, resistance) %>%
  mutate(resistance = as.factor(resistance))

genes <- qpcr %>% select(-resistance) %>% colnames()


lapply(genes, function(gen) {
  
  gp <- ggplot(qpcr, aes(x = resistance, y = get(gen),
                         color = resistance)) +
    geom_violin(size = 2, show.legend = F) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F,) +
    geom_dotplot(binaxis= "y", 
                 stackdir = "center", dotsize = 1, 
                 fill = 1, color = 1) +
    ylab("qPCR expression") +
    xlab("") +
    ggtitle(paste0(toupper(gen), " expression")) +
    theme_bw(base_size = 24) +
    stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                       label.x.npc = "center") +
    scale_x_discrete(labels=c("Sensitive","Resistant")) +
    scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
  
  png(filename = paste0("plots_qPCR/", toupper(gen), ".png"), width = 800, height = 600)
  print(gp)
  dev.off()
  
  gp <- ggplot(qpcr, aes(x = resistance, y = get(gen),
                         color = resistance)) +
    geom_boxplot(size = 2, show.legend = F) +
    geom_dotplot(binaxis= "y", 
                 stackdir = "center",
                 dotsize = 1,  fill = 1, color = 1) +
    ylab("qPCR expression") +
    xlab("") +
    ggtitle(paste0(toupper(gen), " expression")) +
    theme_bw(base_size = 24) +
    stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                       label.x.npc = "center") +
    scale_x_discrete(labels=c("Sensitive","Resistant")) +
    scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
  
  png(filename = paste0("plots_qPCR/", toupper(gen), "_boxplot.png"), width = 800, height = 600)
  print(gp)
  dev.off()
})


combns <- expand_grid(gen1 = genes, gen2 = genes) %>% 
  filter(gen1 != gen2) %>% 
  mutate(name = paste0(pmin(gen1, gen2), "-", pmax(gen1, gen2))) %>%
  distinct(name, .keep_all = T) %>% select (-name)

all_dfs <- mapply(function(gen1, gen2) {
  
  df_qpcr <- qpcr %>% select(expr = all_of(gen1), resistance) %>%
    mutate(gen = toupper(gen1)) %>%
    bind_rows(qpcr %>% select(expr = all_of(gen2), resistance) %>%
                mutate(gen = toupper(gen2)))

  gp <- ggplot(df_qpcr, aes(x = resistance, y = expr)) +
    geom_violin(color="grey", size = 2, show.legend = F) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F,) +
    geom_dotplot(aes(fill = gen), binaxis= "y",
                 stackdir = "center", dotsize = 1) +
    ylab("qPCR expression") +
    xlab("") +
    ggtitle(paste0(toupper(gen1), "-", toupper(gen2), " expression")) +
    theme_bw(base_size = 30) +
    stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                       label.x.npc = "center") +
    scale_x_discrete(labels=c("Sensitive","Resistant")) +
    scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant")) +
    scale_fill_brewer(palette = "Pastel1", name = "")
  
  png(filename = paste0("plots_qPCR/", toupper(gen1), "_", toupper(gen2),  ".png"), 
      width = 800, height = 600)
  print(gp)
  dev.off()
  
  gp <- ggplot(df_qpcr, aes(x = resistance, y = expr)) +
    geom_boxplot(color="grey", size = 2, show.legend = F) +
    geom_dotplot(aes(fill = gen), binaxis= "y",
                 stackdir = "center", dotsize = 1) +
    ylab("qPCR expression") +
    xlab("") +
    ggtitle(paste0(toupper(gen1), "-", toupper(gen2), " expression")) +
    theme_bw(base_size = 30) +
    stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                       label.x.npc = "center") +
    scale_x_discrete(labels=c("Sensitive","Resistant")) +
    scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant")) +
    scale_fill_brewer(palette = "Pastel1", name = "")
  
  png(filename = paste0("plots_qPCR/", toupper(gen1), "_", toupper(gen2),  "_boxplot.png"), 
      width = 800, height = 600)
  print(gp)
  dev.off()
}, combns$gen1, combns$gen2, SIMPLIFY = FALSE)
