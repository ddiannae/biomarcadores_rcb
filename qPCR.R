library(readr)
library(dplyr)
#library(tidyr)
library(ggplot2)
library(ggpubr)
library(janitor)

qpcr <- read_csv("data_qPCR/QPCR.csv") %>%
  clean_names() %>%
  mutate(resistant = as.factor(resistant))

dir.create("plots_qPCR")

gp <- ggplot(qpcr, aes(x = resistant, y = ndufaf3,
                           color = resistant)) +
  geom_boxplot() +
  #stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F,) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("qPCR expression") +
  xlab("") +
  ggtitle("NDUFAF3 expression") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))
  
png(filename = "plots_qPCR/NDUFAF3_boxplot.png", width = 800, height = 400)
print(gp)
dev.off()

gp <- ggplot(qpcr, aes(x = resistant, y = gria4,
                       color = resistant)) +
  geom_boxplot() +
  #stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("qPCR expression") +
  xlab("") +
  ggtitle("GRIA4 expression") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank())  +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots_qPCR/GRIA4_boxplot.png", width = 800, height = 400)
print(gp)
dev.off()

gp <- ggplot(qpcr, aes(x = resistant, y = slc12a,
                       color = resistant)) +
  geom_boxplot() +
  #stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("qPCR expression") +
  xlab("") +
  ggtitle("SLC12A expression") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "top",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank())  +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots_qPCR/SLC12A_boxplot.png", width = 800, height = 400)
print(gp)
dev.off()
