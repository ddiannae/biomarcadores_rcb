library(readr)
library(dplyr)
library(janitor)
library(purrr)
library(broom)

qpcr <- read_csv("data_qPCR/originales_qPCR.csv") %>%
  clean_names() 

datos_mama <- read_csv("data_qPCR/base_completa.csv") %>%
  clean_names() 

qpcr <- qpcr %>% inner_join(datos_mama, by = "folio") %>%
  filter(q_pcr32 == 1)
## Son las variables que seleccionaste
qpcr <- qpcr %>% select(recurrence, death, resistant, age, c_t, c_n, tnm, er,
                her2, ndufaf3, gria4, slc12a)

qpcr <- qpcr %>% 
  mutate(recurrence = as.factor(recurrence),
         death = as.factor(death),
         resistant = as.factor(resistant),
         tnm = as.numeric(grepl("^III", tnm)),
         ndufaf3 = if_else(ndufaf3 < median(ndufaf3, na.rm = TRUE), 0, 1),
         gria4 = if_else(gria4 < median(gria4, na.rm = TRUE), 0, 1),
         slc12a = if_else(slc12a < median(slc12a, na.rm = TRUE), 0, 1))


### Modelos univariados 

## variables independientes
ivars <- c("recurrence", "death", "resistant")
## variables dependientes
dvars <- colnames(qpcr)[!colnames(qpcr) %in% ivars]

recurrence_stats <- map_chr(dvars, ~paste0(ivars[1], "~", .x)) %>%
    map(.f = ~glm(formula = as.formula(.x),      
    family = "binomial", data = qpcr, na.action = na.omit)) %>% 
    map_df(.f = ~tidy(.x,  conf.int = TRUE)) %>%
    filter(term != "(Intercept)")

recurrence_stats %>% filter(p.value < 0.05)
# A tibble: 0 × 7

death_stats <- map_chr(dvars, ~paste0(ivars[2], "~", .x)) %>%
  map(.f = ~glm(formula = as.formula(.x),      
                family = "binomial", data = qpcr, na.action = na.omit)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE)) %>%
  filter(term != "(Intercept)")

death_stats %>% filter(p.value < 0.05)
# # A tibble: 1 × 7
# term  estimate std.error statistic p.value conf.low conf.high
# <chr>    <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 age      0.206    0.0996      2.07  0.0388   0.0420     0.451


resistant_stats <- map_chr(dvars, ~paste0(ivars[3], "~", .x)) %>%
  map(.f = ~glm(formula = as.formula(.x),      
                family = "binomial", data = qpcr, na.action = na.omit)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE)) %>%
  filter(term != "(Intercept)")

resistant_stats %>% filter(p.value < 0.05)
# A tibble: 1 × 7
# term   estimate std.error statistic p.value conf.low conf.high
# <chr>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 slc12a    -1.70     0.801     -2.13  0.0332    -3.38    -0.199

### Modelos multivariados
data_resistant <- qpcr %>% select(-recurrence, -death) 
multi_resistant_stats <-  glm(resistant ~., data = data_resistant,
                              family="binomial") %>% tidy(conf.int = TRUE)
multi_resistant_stats %>% filter(p.value < 0.05)
# A tibble: 0 × 7

data_death <- qpcr %>% select(-recurrence, -resistant) 
multi_death_stats <-  glm(death ~., data = data_death,
                              family="binomial") %>% tidy(conf.int = TRUE)
multi_death_stats %>% filter(p.value < 0.05)
# A tibble: 0 × 7

data_recurrence <- qpcr %>% select(-death, -resistant) 
multi_recurrence_stats <-  glm(recurrence ~., data = data_recurrence,
                          family="binomial") %>% tidy(conf.int = TRUE)
multi_recurrence_stats %>% filter(p.value < 0.05)
# A tibble: 0 × 7
