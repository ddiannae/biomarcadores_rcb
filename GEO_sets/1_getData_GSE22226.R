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
dir.create("data_GSE22226")

# Obtenemos los datos del set de GEO
gse <- getGEO("GSE22226", GSEMatrix=FALSE)

# Analizando la primera muestra. Con Meta se obtiene la información asociada a la muestra
sample1 <- Meta(GSMList(gse)[[1]])
str(sample1)

gse_meta <- lapply(GSMList(gse), Meta)

# Obtenemos la información relevante para todas las muestras
# characteristics_ch2 es la cuarta entrada en la lista
gsms <- gse_meta %>% map(4)

# Convertimos la información de cada muestra a un data frame
gsms <- lapply(seq_along(gsms), function(i) {
  gsm <- gsms[[i]]
  # Separamos los datos en los campos, por ejemplo 'sample id: 1002'
  all_chras <- lapply(str_split(gsm, ":"), str_trim)
  
  tchrs <- lapply(all_chras, function(char) {
    ifelse(length(char) == 3, 
           return(list(attr = paste(str_replace_all(char[[1]], " ", "_"),
                                    str_replace_all(char[[2]], " ", "_"), "_"),
                       value = char[[3]])),
           return(list(attr = str_replace_all(char[[1]], " ", "_"),
                       value = char[[2]])))
  })
  
  # Construimos un data frame con los pares de dato y valor que obtuvimos al separar por :
  tchrs <- bind_rows(tchrs)
  
  # Los datos que faltan dicen 'NA', asi que los convertimos a NA
  tchrs[tchrs == "NA"] <- NA
  
  # Le ponemos el mismo id a todo los campos asociados a la muestra
  tchrs$id <- i
  return(tchrs)
})

# Unimos los datos de las muestras y convertimos el data frame a formato ancho
gsms <- bind_rows(gsms) %>%
  pivot_wider(id_cols = id, names_from = attr, values_from = value) %>%
  clean_names() %>%
  select(-id, -tissue)

# Las primeras entradas del data frame
head(gsms)

# Creamos otro data frame con los nombres de las muestras en GEO, como GSM615096
gsms_names <- gse_meta %>%
  map(`[`, c("geo_accession", "description", "platform_id")) %>% 
  bind_rows()

# Lo unimos al dataframe con todos los datos
gsms <- gsms %>% 
  inner_join(gsms_names, by = c("experiment_name" = "description")) %>%
  select(geo_accession, platform_id, experiment_name, everything())

tibble(column = colnames(gsms)) %>%
  write_tsv("data_GSE22226/columns_description.tsv")

colnames(gsms) = c(
  "geo_accession",
  "platform_id",
  "experiment_name",
  "slide_name",
  "i_spy_id",
  "study",
  "age",
  "histologic_grade",
  "histology",
  "clinical_tumor_size",
  "clinical_t_stage",
  "er_status",
  "pgr_status",
  "pam50_subtype",
  "neoadjuvant_chemotherapy",
  "pathological_complete_response_pcr",
  "rcb_class",
  "rfs_time",
  "rfs_indicator",
  "os_time",
  "survival_status",
  "her2_status")

# Seleccionamos las muestras de pacientes localmente avanzadas                 
local_avan <- gsms %>% 
  filter(er_status == 1 & clinical_t_stage %in% 2:3)

# Agregamos el campo rs_group, para asignar R de resistente o S de sensible, según lo 
# indicado en el campo rcs_class
local_avan <- local_avan %>% 
  mutate(rs_group = ifelse(is.na(rcb_class) |rcb_class == "Unavailable" , NA, 
                              ifelse(str_detect(rcb_class, "II|II"), "R", "S"))) %>%
  select(geo_accession, platform_id, experiment_name, rs_group, everything())


# Las primeras entradas del data frame
head(local_avan)

# Guardamos un archivo con datos de las pacientes resistentes
local_avan %>% filter(rs_group =="R") %>% write_tsv("data_GSE22226/resistentes.tsv")

# Guardamos un archivo con datos de las pacientes sensibles
local_avan %>% filter(rs_group =="S") %>% write_tsv("data_GSE22226/sensibles.tsv")

# Guardamos el data frame con todos los datos
local_avan %>% write_tsv("data_GSE22226/localmente_avanzadas.tsv")

# El estudio tiene dos GPL, que es el GEO con la descripción del arreglo
gse@gpls[[1]]@header$geo_accession
gse@gpls[[2]]@header$geo_accession

GPL1708 <- gse@gpls[[1]]@dataTable@table %>% clean_names()  %>% 
  select(id, col, row, name, spot_id, control_type, gene_symbol) %>%
  filter(gene_symbol != "")

GPL4133 <- gse@gpls[[2]]@dataTable@table %>% clean_names()  %>% 
  select(id, col, row, name, spot_id_5, control_type, gene_symbol) %>%
  filter(gene_symbol != "")

## Buscamos los genes de ambas plataformas
plat_inter <- intersect(GPL1708$gene_symbol, GPL4133$gene_symbol)

GPL1708 <- GPL1708 %>% filter(gene_symbol %in% plat_inter)
GPL4133 <- GPL4133 %>% filter(gene_symbol %in% plat_inter)

## Guardamos ambas listas de features
GPL1708 %>% write_tsv("data_GSE22226/GPL1708_features.tsv")
GPL4133 %>% write_tsv("data_GSE22226/GPL4133_features.tsv")

local_avan_gsms <- local_avan %>% filter(!is.na(rs_group)) %>% 
  pull(geo_accession)
## Ahora obtenemos los datos de expresión
gsms_data <- lapply(gse@gsms[local_avan_gsms], function(gsm) {
  dt <- GEOquery::dataTable(gsm)@table %>% clean_names() %>%
    select(id_ref, value)
  dt$geo_accession <- GEOquery::Meta(gsm)$geo_accession
  return(dt)
})
gsms_data <- bind_rows(gsms_data) %>% 
  pivot_wider(id_cols = id_ref, names_from = geo_accession, 
              values_from = value)

gsms_data_GPL1708 <- gsms_data %>% select(id_ref, 
  any_of(local_avan %>% filter(platform_id == "GPL1708") %>% 
           pull(geo_accession))) %>%
  filter(id_ref %in% (GPL1708 %>% pull(id)))

gsms_data_GPL4133 <- gsms_data %>% select(id_ref, 
  any_of(local_avan %>% filter(platform_id == "GPL4133") %>% 
           pull(geo_accession))) %>%
  filter(id_ref %in% (GPL4133 %>% pull(id)))

gsms_data_GPL1708 %>% write_tsv("data_GSE22226/GPL1708_expression.tsv")
gsms_data_GPL4133 %>% write_tsv("data_GSE22226/GPL4133_expression.tsv")
