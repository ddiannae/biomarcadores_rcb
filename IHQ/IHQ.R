library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(janitor)

ihq <- read_excel("../data/Datos crudos.xlsx", sheet = "IHQ_n=31ER+", skip=3, 
          col_names = c("folio", "resistant", "NDUFAF3_porcentaje_positividad",
                        "NDUFAF3_positividad", "NDUFAF3_intensidad", "NDUFAF3_score", 
                        "GRIA4_porcentaje_positividad",  "GRIA4_positividad", 
                        "GRIA4_intensidad", "GRIA4_score")) 


### NDUFAF3 violin y boxplots
png(filename = "plots_IHQ/NDUFAF3_score_violin.png", width = 800, height = 600)
ihq %>% select(resistant, NDUFAF3_score) %>% 
  filter(!is.na(NDUFAF3_score)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = NDUFAF3_score, color = resistant)) +
  geom_violin(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 1,  fill = 1, color = 1) +
  stat_summary(fun = "median", geom = "crossbar", show.legend = F) +
  ylab("IHC Score") +
  xlab("") +
  ggtitle("NDUFAF3 score") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/NDUFAF3_score_boxplot.png", width = 800, height = 600)
ihq %>% select(resistant, NDUFAF3_score) %>% 
  filter(!is.na(NDUFAF3_score)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = NDUFAF3_score, color = resistant)) +
  geom_boxplot(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 1,  fill = 1, color = 1) +
  ylab("IHC Score") +
  xlab("") +
  ggtitle("NDUFAF3 score") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/NDUFAF3_positividad_violin.png", width = 800, height = 600)
ihq %>% select(resistant, NDUFAF3_positividad) %>% 
  filter(!is.na(NDUFAF3_positividad)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = NDUFAF3_positividad, color = resistant)) +
  geom_violin(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 0.4,  fill = 1, color = 1) +
  stat_summary(fun = "median", geom = "crossbar", show.legend = F) +
  ylab("IHC positivity score") +
  xlab("") +
  ggtitle("NDUFAF3 positivity") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/NDUFAF3_positividad_boxplot.png", width = 800, height = 600)
ihq %>% select(resistant, NDUFAF3_positividad) %>% 
  filter(!is.na(NDUFAF3_positividad)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = NDUFAF3_positividad, color = resistant)) +
  geom_boxplot(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 0.4,  fill = 1, color = 1) +
  ylab("IHC positivity score") +
  xlab("") +
  ggtitle("NDUFAF3 positivity") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()


png(filename = "plots_IHQ/NDUFAF3_intensity_violin.png", width = 800, height = 600)
ihq %>% select(resistant, NDUFAF3_intensidad) %>% 
  filter(!is.na(NDUFAF3_intensidad)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = NDUFAF3_intensidad, color = resistant)) +
  geom_violin(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 0.3,  fill = 1, color = 1) +
  stat_summary(fun = "median", geom = "crossbar", show.legend = F) +
  ylab("IHC intensity score") +
  xlab("") +
  ggtitle("NDUFAF3 intensity") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/NDUFAF3_intensity_boxplot.png", width = 800, height = 600)
ihq %>% select(resistant, NDUFAF3_intensidad) %>% 
  filter(!is.na(NDUFAF3_intensidad)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = NDUFAF3_intensidad, color = resistant)) +
  geom_boxplot(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 0.3,  fill = 1, color = 1) +
  ylab("IHC intensity score") +
  xlab("") +
  ggtitle("NDUFAF3 intensity") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

#### GRIA4 violin y boxplots
png(filename = "plots_IHQ/GRIA4_score_violin.png", width = 800, height = 600)
ihq %>% select(resistant, GRIA4_score) %>% 
  filter(!is.na(GRIA4_score)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = GRIA4_score, color = resistant)) +
  geom_violin(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 1,  fill = 1, color = 1) +
  stat_summary(fun = "median", geom = "crossbar", show.legend = F) +
  geom_signif(comparisons = list(c("0", "1")), 
              map_signif_level=TRUE, size = 1, color = "black", 
              vjust = 0.5, textsize = 10) +
  ylab("IHC Score") +
  xlab("") +
  ggtitle("GRIA4 score") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/GRIA4_score_boxplot.png", width = 800, height = 600)
ihq %>% select(resistant, GRIA4_score) %>% 
  filter(!is.na(GRIA4_score)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = GRIA4_score, color = resistant)) +
  geom_boxplot(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 1,  fill = 1, color = 1) +
  geom_signif(comparisons = list(c("0", "1")), 
              map_signif_level=TRUE, size = 1, color = "black", 
              vjust = 0.5, textsize = 10) +
  ylab("IHC Score") +
  xlab("") +
  ggtitle("NDUFAF3 score") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/GRIA4_positividad_violin.png", width = 800, height = 600)
ihq %>% select(resistant, GRIA4_positividad) %>% 
  filter(!is.na(GRIA4_positividad)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = GRIA4_positividad, color = resistant)) +
  geom_violin(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 0.4,  fill = 1, color = 1) +
  stat_summary(fun = "median", geom = "crossbar", show.legend = F) +
  ylab("IHC positivity score") +
  xlab("") +
  ggtitle("GRIA4 positivity") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  geom_signif(comparisons = list(c("0", "1")), 
              map_signif_level=TRUE, size = 1, color = "black", 
              vjust = 0.5, textsize = 10) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/GRIA4_positividad_boxplot.png", width = 800, height = 600)
ihq %>% select(resistant, GRIA4_positividad) %>% 
  filter(!is.na(GRIA4_positividad)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = GRIA4_positividad, color = resistant)) +
  geom_boxplot(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 0.4,  fill = 1, color = 1) +
  ylab("IHC positivity score") +
  xlab("") +
  ggtitle("GRIA4 positivity") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  geom_signif(comparisons = list(c("0", "1")), 
              map_signif_level=TRUE, size = 1, color = "black", 
              vjust = 0.5, textsize = 10) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/GRIA4_intensidad_violin.png", width = 800, height = 600)
ihq %>% select(resistant, GRIA4_intensidad) %>% 
  filter(!is.na(GRIA4_intensidad)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = GRIA4_intensidad, color = resistant)) +
  geom_violin(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 0.3,  fill = 1, color = 1) +
  stat_summary(fun = "median", geom = "crossbar", show.legend = F) +
  ylab("IHC intensity score") +
  xlab("") +
  ggtitle("GRIA4 intensity") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()

png(filename = "plots_IHQ/GRIA4_intensidad_boxplot.png", width = 800, height = 600)
ihq %>% select(resistant, GRIA4_intensidad) %>% 
  filter(!is.na(GRIA4_intensidad)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  ggplot(., aes(x = resistant, y = GRIA4_intensidad, color = resistant)) +
  geom_boxplot(size = 2, show.legend = F) +
  geom_dotplot(binaxis= "y", binwidth = 1/2,
               stackdir = "center",
               dotsize = 0.3,  fill = 1, color = 1) +
  ylab("IHC intensity score") +
  xlab("") +
  ggtitle("GRIA4 intensity") +
  theme_bw(base_size = 30) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center", ) +
  scale_x_discrete(labels=c("Sensitive","Resistant")) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
dev.off()