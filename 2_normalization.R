# Para instalar los paquetes
# install.packages(c("stringr", "readr", "dplyr", "purrr", "reshape2", "ggplot2"))
library(stringr)
library(readr)
library(dplyr)
library(purrr)
library(reshape2)
library(ggplot2)

# Para instalar paquetes de bioconductor
# BiocManager::install("affyPLM")
# BiocManager::install("gcrma")
library(affyPLM)
library(gcrma)

# Creamos el directorio para guardar las gráficas
dir.create("plots")

# Leemos el data frame con la información de las muestras
local_avan <- read_tsv("data/localmente_avanzadas.tsv")

# local_avan <- read_tsv("data/localmente_avanzadas_pCR_RD.tsv")

# Listamos los CEL files descargados en el paso anterior
cel_files <- list.files(path = "raw", pattern = "*/*.CEL.gz", recursive = TRUE, full.names = T)

# cel_files <- list.files(path = "raw_pcr_rd", pattern = "*/*.CEL.gz", recursive = TRUE, full.names = T)

# Leemos los cel files en un objeto affybatch
cels <- ReadAffy(filenames = cel_files, 
                 # Los nombres de las muestras se obtienen del nombre del archivo
                 sampleNames = unlist(map(strsplit(unlist(map(strsplit(cel_files, split = "\\/"), 3)), "_"), 1)))

# Filtramos el data frame para que solo sean las muestras con CEL files
# Se quitan las muestras que no tenían información de RCB-class
local_avan <- local_avan %>% filter(name %in% affy::sampleNames(cels))

# Función para guardar boxplots
saveBoxPlot <- function(data, name) {
  
  # Convertimos el nombre de la gráfica a lower case y unido con _ para el nombre del archivo
  low_name <- tolower(str_replace_all(name, " ", "_"))
  
  # Guardamos la gráfica
  png(filename = paste0("plots/", low_name, "_boxplot.png"), width = 800, height = 400)
  boxplot(data, names = NA, main = name)
  dev.off()
}

# Función para guardar density plots
saveDensityPlot <- function(data, name, group_data) {
  
  # Convertimos el nombre de la gráfica a lower case y unido con _ para el nombre del archivo
  low_name <- tolower(str_replace_all(name, " ", "_"))
  
  # Convertimos los datos a formato largo y unimos el data frame con los grupos
  melted_data <- melt(data)
  melted_data <- melted_data %>% inner_join(group_data,
                                        by = c("Var2" = "name"))
  
  # Construimos la gráfica
  pl <- ggplot(data = melted_data, aes(x=value, group=Var2, colour=rs_group)) + 
    geom_density() + 
    xlab("log2(expr)") + 
    ggtitle(name)
  
  # Guardamos la gráfica
  png(filename = paste0("plots/", low_name, "_density.png"), width = 800, height = 400)
  print(pl)
  dev.off()
  
}

# Guardamos gráficas de los datos de expresión crudos
saveBoxPlot(cels, "Raw expression")
saveDensityPlot(log2(exprs(cels)), "Raw expression", local_avan %>% 
                  dplyr::select(name, rs_group))

# saveBoxPlot(cels, "Raw expression pCR RD")
# saveDensityPlot(log2(exprs(cels)), "Raw expression pCR RD", local_avan %>% 
#                   dplyr::select(name, rs_group))

# Normalización por RMA
rma_data <- rma(cels)
saveBoxPlot(rma_data, "RMA normalized")
# Los resultados de RMA ya están en log2
saveDensityPlot(exprs(rma_data), "RMA normalized", local_avan %>% 
                  dplyr::select(name, rs_group))

saveBoxPlot(rma_data, "RMA normalized pCR RD")
saveDensityPlot(exprs(rma_data), "RMA normalized pCR RD", local_avan %>% 
                  dplyr::select(name, rs_group))

# Guardamos los datos de expresión, con todo y el nombre de los features
expr <- as.data.frame(exprs(rma_data))
expr$feature <- rownames(expr)
expr %>% dplyr::select(feature, everything())%>% 
  write_tsv("data/rs_rma_normalized.tsv")

# expr %>% dplyr::select(feature, everything())%>% 
#   write_tsv("data/rs_rma_normalized_pcr_rd.tsv")

# Normalización por GCRMA
gcrma_data <- gcrma(cels)
saveBoxPlot(gcrma_data, "GCRMA normalized ")
# Los resultados de RMA ya están en log2
saveDensityPlot(exprs(gcrma_data), "GCRMA normalized", local_avan %>% 
                  dplyr::select(name, rs_group))

# saveBoxPlot(gcrma_data, "GCRMA normalized pCR RD")
# # Los resultados de RMA ya están en log2
# saveDensityPlot(exprs(gcrma_data), "GCRMA normalized pCR RD", local_avan %>% 
#                   dplyr::select(name, rs_group))

# Guardamos los datos de expresión,  con todo y el nombre de los features
expr <- as.data.frame(exprs(gcrma_data))
expr$feature <- rownames(expr)
expr %>% dplyr::select(feature, everything())%>% 
  write_tsv("data/rs_gcrma_normalized.tsv")

# expr %>% dplyr::select(feature, everything())%>% 
#   write_tsv("data/rs_gcrma_normalized_pcr_rd.tsv")
