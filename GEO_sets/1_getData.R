# Para instalar los paquetes
# install.packages(c("stringr", "dplyr", "purrr", "readr", "tidyr"))
library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(janitor)

# Para instalar los paquetes de Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")
library(GEOquery)

# Creamos directorio para guardar los datos que obtenemos
dir.create("data_la_extra")

# Obtenemos los datos del set de GEO
gse <- getGEO("GSE25066", GSEMatrix=FALSE)

# El estudio solo tiene un GPL, que es el GEO con la descripción del arreglo
# Guardamos los datos de las features
gse@gpls[[1]]@dataTable@table %>% clean_names() %>%
  select(id, gb_acc, spot_id, gene_title, gene_symbol, entrez_gene_id) %>%
  write_tsv("data_la_extra/features.tsv")

# Analizando la primera muestra. Con Meta se obtiene la información asociada a la muestra
sample1 <- Meta(GSMList(gse)[[1]])
str(sample1)

# Obtenemos la información relevante para todas las muestras
# characteristics_ch1 es la segunda entrada en la lista
gsms <- lapply(GSMList(gse), Meta) %>% map(2)

# Convertimos la información de cada muestra a un data frame
gsms <- lapply(gsms, function(gsm) {
  
  # Separamos los datos en los campos, por ejemplo 'sample id: 1002'
  all_chrs <- lapply(str_split(gsm, ":"), str_trim)
  
  # Construimos un data frame con los pares de dato y valor que obtuvimos al separar por :
  tchrs <- tibble(attr = str_replace_all(map_chr(all_chrs, 1), " ", "_"), value = map_chr(all_chrs, 2))
  
  # Los datos que faltan dicen 'NA', asi que los convertimos a NA
  tchrs[tchrs == "NA"] <- NA
  
  # Le ponemos el mismo id a todo los campos asociados a la muestra
  tchrs$id <- tchrs %>% filter(attr == "sample_id") %>% pull(value)
  return(tchrs)
})

# Unimos los datos de las muestras y convertimos el data frame a formato ancho
gsms <- bind_rows(gsms) %>%
  pivot_wider(id_cols = id, names_from = attr, values_from = value) %>%
  select(-id)

# Las primeras entradas del data frame
head(gsms)

# Creamos otro data frame con los nombres de las muestras en GEO, como GSM615096
gsms_names <- tibble(name = names(GSMList(gse)), 
                     sample_id = lapply(GSMList(gse), Meta) %>% 
                       map("title") %>% unlist())

# Lo unimos al dataframe con todos los datos
gsms <- gsms %>% inner_join(gsms_names, by ="sample_id") %>%
  select(name, sample_id, everything())

# Seleccionamos las muestras de pacientes localmente avanzadas                 
local_avan <- gsms %>% 
  filter(er_status_ihc == "P" & 
           clinical_ajcc_stage %in% c("IIA", "IIIA", "IIB", "IIIB", "IIIC"))

# Agregamos el campo rs_group, para asignar R de resistente o S de sensible, según lo 
# indicado en el campo rcs_class
local_avan <- local_avan %>% 
  mutate(rs_group = ifelse(is.na(pathologic_response_rcb_class), NA, 
                              ifelse(pathologic_response_rcb_class %in% c("RCB-II", "RCB-III"), "R", "S"))) %>%
  select(name, sample_id, rs_group, everything()) %>% 
  filter(!is.na(rs_group))


# Para filtrar por PCR RD
# local_avan <- local_avan %>% 
#   mutate(rs_group = ifelse(is.na(pathologic_response_pcr_rd), NA, 
#                            ifelse(pathologic_response_pcr_rd == "RD", "R", "S"))) %>%
#   select(name, sample_id, rs_group, everything())

# Las primeras entradas del data frame
head(local_avan)

# Guardamos un archivo con datos de las pacientes resistentes
local_avan %>% filter(rs_group =="R") %>% write_tsv("data_la_extra/resistentes.tsv")

# Para almacenar RD
# local_avan %>% filter(rs_group =="R") %>% write_tsv("data/resistentes_RD.tsv")

# Guardamos un archivo con datos de las pacientes sensibles
local_avan %>% filter(rs_group =="S") %>% write_tsv("data_la_extra/sensibles.tsv")

# Para almacenar pCR
# local_avan %>% filter(rs_group =="S") %>% write_tsv("data/sensibles_pCR.tsv")

# Guardamos el data frame con todos los datos
local_avan %>% write_tsv("data_la_extra/localmente_avanzadas.tsv")

# local_avan %>% write_tsv("data/localmente_avanzadas_pCR_RD.tsv")

# Creamos directorio para guardar los CEL files
dir.create("raw_la_extra")

# dir.create("raw_pcr_rd")

# Descargamos los CEL files usando los nombres en el data frame
local_avan %>% filter(rs_group %in% c("R", "S")) %>% 
  pull(name) %>% map(getGEOSuppFiles, baseDir = "raw_la_extra")

# local_avan %>% filter(rs_group %in% c("R", "S")) %>% 
#   pull(name) %>% map(getGEOSuppFiles, baseDir = "raw_pcr_rd")
