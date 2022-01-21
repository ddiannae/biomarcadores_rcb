library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(janitor)

ihq <- read_csv("data_IHQ/IHQ.csv") %>%
  clean_names()

# a) Comparación de la expresión entre sensibles y resistentes
# antes de recibir la quimio: IHC_score_NDUFpreNAC (sensibles vs resistentes)
pre_nduf <- ihq %>% select(ihc_score_ndu_fpre_nac, score_positivity_ndu_fpre_nac, 
                           score_intensity_ndu_fpre_nac,
                           rcb_class, resistant) %>%
  filter(!is.na(ihc_score_ndu_fpre_nac)) %>%
  mutate(resistant = as.factor(resistant))

gp <- ggplot(pre_nduf, aes(x = resistant, y = ihc_score_ndu_fpre_nac,
                           color = resistant)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("IHC Score") +
  xlab("") +
  ggtitle("NDUFAF3 score preNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F,  label.y.npc = "bottom",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots_IHQ/NDUFAF3_score.png", width = 800, height = 400)
print(gp)
dev.off()

gp <- ggplot(pre_nduf, aes(x = resistant, y = score_positivity_ndu_fpre_nac,
                           color = resistant)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("IHC Positivity") +
  xlab("") +
  ggtitle("NDUFAF3 score preNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F, label.y.npc = "bottom",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots_IHQ/NDUFAF3_positivity.png", width = 800, height = 400)
print(gp)
dev.off()

gp <- ggplot(pre_nduf, aes(x = resistant, y = score_intensity_ndu_fpre_nac,
                           color = resistant)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("IHC Intensity") +
  xlab("") +
  ggtitle("NDUFAF3 score preNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F, label.y.npc = "bottom",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots_IHQ/NDUFAF3_intensity.png", width = 800, height = 400)
print(gp)
dev.off()

# b) Comparación de la expresión entre resistentes antes y después de 
# recibir la quimio: IHC_score_NDUFpreNAC vs  IHC_score_NDUFposNAC 
# (seleccionando solo a resistentes)
r_nduf <- ihq %>% select(ic_bx1, ihc_score_ndu_fpre_nac, 
                         score_positivity_ndu_fpre_nac, 
                         score_intensity_ndu_fpre_nac, rcb_class, resistant) %>%
      filter(resistant == 1 & !is.na(ihc_score_ndu_fpre_nac)) %>% 
      rename("score" = "ihc_score_ndu_fpre_nac",
             "positivity" = "score_positivity_ndu_fpre_nac",
             "intensity" = "score_intensity_ndu_fpre_nac") %>%
      mutate(class = "preNac") %>%
      bind_rows(ihq %>% select(ic_bx1, ihc_score_ndu_fpos_nac, 
                         score_positivity_ndu_fpos_nac, 
                         score_intensity_ndu_fpos_nac, rcb_class, resistant) %>%
                   filter(resistant == 1 & !is.na(ihc_score_ndu_fpos_nac)) %>% 
                   rename("score" = "ihc_score_ndu_fpos_nac",
                          "positivity" = "score_positivity_ndu_fpos_nac",
                          "intensity" = "score_intensity_ndu_fpos_nac")  %>%
                  mutate(class = "posNac")) %>%
      mutate(class = factor(class, levels = c("preNac", "posNac")))
  

gp <- ggplot(r_nduf, aes(x = class, y = score, color = class)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  
  labs(x = "", y = "IHC Score", 
       title = "Resistant NDUFAF3 score preNAC vs posNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F, label.y.npc = "bottom",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set1", name = "")

png(filename = "plots_IHQ/resistant_NDUFAF3_score.png", width = 800, height = 400)
print(gp)
dev.off()

gp <- ggplot(r_nduf, aes(x = class, y = intensity, color = class)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  labs(x = "", y = "IHC Intensity", 
       title = "Resistant NDUFAF3 intensity preNAC vs posNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F, label.y.npc = "bottom",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set1", name = "")

png(filename = "plots_IHQ/resistant_NDUFAF3_intensity.png", width = 800, height = 400)
print(gp)
dev.off()

gp <- ggplot(r_nduf, aes(x = class, y = positivity, color = class)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  labs(x = "", y = "IHC Positivity", 
       title = "Resistant NDUFAF3 positivity preNAC vs posNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F, label.y.npc = "bottom",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set1", name = "")

png(filename = "plots_IHQ/resistant_NDUFAF3_positivity.png", width = 800, height = 400)
print(gp)
dev.off()

# a) Comparación de la expresión entre sensibles y resistentes
# antes de recibir la quimio: IHC_score_NDUFpreNAC (sensibles vs resistentes)
pre_gria <- ihq %>% filter(pre_nac == 1) %>%
  select(folio, ihc_score_gri_apre_nac, score_positivity_gri_apre_nac, 
                           score_intensity_gri_apre_nac,
                           rcb_class, resistant) %>%
  filter(!is.na(ihc_score_gri_apre_nac)) %>%
  mutate(resistant = as.factor(resistant)) %>%
  arrange(resistant, ihc_score_gri_apre_nac)


gp <- ggplot(pre_gria, aes(x = resistant, y = ihc_score_gri_apre_nac,
                           color = resistant)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("IHC Score") +
  xlab("") +
  ggtitle("GRIA4 score preNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F, label.y.npc = "bottom",  
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots_IHQ/GRIA4_score.png", width = 800, height = 400)
print(gp)
dev.off()

gp <- ggplot(pre_gria, aes(x = resistant, y = score_intensity_gri_apre_nac,
                           color = resistant)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("IHC Positivity") +
  xlab("") +
  ggtitle("GRIA4 positivity preNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F, label.y.npc = "bottom", 
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots_IHQ/GRIA4_positivity.png", width = 800, height = 400)
print(gp)
dev.off()


gp <- ggplot(pre_gria, aes(x = resistant, y = score_intensity_gri_apre_nac,
                           color = resistant)) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, show.legend = F) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  ylab("IHC Intensity") +
  xlab("") +
  ggtitle("GRIA4 intensity preNAC") +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test", show.legend = F, label.y.npc = "bottom",
                     label.x.npc = "center") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set2", name = "", labels = c("Sensitive", "Resistant"))

png(filename = "plots_IHQ/GRIA4_intensity.png", width = 800, height = 400)
print(gp)
dev.off()

