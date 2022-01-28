library(readr)
library(dplyr)

deg_lnc <- read_tsv("data_lncrna/deg.tsv")
deg_geo <- read_tsv("data/deg_sensibles_vs_resistentes.tsv")
load("data_lncrna/deg_originales.RData")

deg_geo <- deg_geo %>% filter(adj_p_val < 0.05)

brca_filtered <- brca.res.df %>% 
  filter(padj < 0.05) %>% as_tibble()

gene_intersect <- brca_filtered %>% inner_join(deg_geo, by = c("hgnc_symbol" = "gene_symbol" ))

deg_geo_extra <- read_tsv("data_la_extra/deg_sensibles_vs_resistentes.tsv")
deg_geo_extra <- deg_geo_extra %>% filter(adj_p_val < 0.05)
gene_intersect_extra <- brca_filtered %>% inner_join(deg_geo_extra, by = c("hgnc_symbol" = "gene_symbol" ))
gene_intersect_extra %>% select(hgnc_symbol, log_fc, p_value, log2FoldChange, pvalue)

  #   hgnc_symbol log_fc   p_value log2FoldChange   pvalue
  #   <chr>        <dbl>     <dbl>          <dbl>    <dbl>
  # 1 SEPHS1       0.387 0.0000883          -1.43 0.000494