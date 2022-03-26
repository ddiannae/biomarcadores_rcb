library(readr)
library(readxl)
library(dplyr)
library(janitor)
library(purrr)
library(broom)

pacientes <- read_excel("../data/Base de datos-Incanet-Patología-INCAN .xlsx", 
                        skip = 2) %>% clean_names()

qpcr <- read_excel("../data/Datos crudos.xlsx", sheet = "qPCR_n=24ER+", 
                   skip = 1) %>% clean_names()

qpcr <- qpcr %>% select(colnames(qpcr)[c(1,35:42)]) 
colnames(qpcr) <- c("folio", "abhd14b", "ndufaf3", "tex264", "ghr", "gria4", 
                    "hgd", "slc12a", "sostcd1")

qpcr <- qpcr %>% inner_join(pacientes, by = "folio") %>%
  mutate(resistance = as.factor(resistance), 
         days = as.numeric(pmin(difftime(last_followup, nac_start, units = "day"),
                                difftime(recurrence_date, nac_start, units = "day"), na.rm = TRUE)))

## Las variables seleccionadas
qpcr <- qpcr %>% select(recurrence, death, resistance, age, c_t, c_n, tnm, luminal_a, 
                luminal_b_her2neg, luminal_b_her2pos,
                ndufaf3, gria4, slc12a)

qpcr <- qpcr %>% 
  mutate(recurrence = as.factor(recurrence),
         death = as.factor(death),
         luminal_a = as.factor(luminal_a),
         luminal_b_her2neg = as.factor(luminal_b_her2neg),
         luminal_b_her2pos = as.factor(luminal_b_her2pos),
         resistance = as.factor(resistance),
         tnm = as.numeric(grepl("^III", tnm)),
         ndufaf3 = if_else(ndufaf3 < median(ndufaf3, na.rm = TRUE), 0, 1),
         gria4 = if_else(gria4 < median(gria4, na.rm = TRUE), 0, 1),
         slc12a = if_else(slc12a < median(slc12a, na.rm = TRUE), 0, 1))


### Modelos univariados 

## variables independientes
ivars <- c("recurrence", "death", "resistance")
## variables dependientes
dvars <- colnames(qpcr)[!colnames(qpcr) %in% ivars]

recurrence_stats <- map_chr(dvars, ~paste0(ivars[1], "~", .x)) %>%
    map(.f = ~glm(formula = as.formula(.x),      
    family = "binomial", data = qpcr, na.action = na.omit)) %>% 
    map_df(.f = ~tidy(.x,  conf.int = TRUE)) %>%
    filter(term != "(Intercept)")

recurrence_stats %>% filter(p.value < 0.05)
# # A tibble: 1 × 7
# term  estimate std.error statistic p.value conf.low conf.high
# <chr>    <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 c_t       2.25      1.08      2.09  0.0367    0.563      5.21

death_stats <- map_chr(dvars, ~paste0(ivars[2], "~", .x)) %>%
  map(.f = ~glm(formula = as.formula(.x),      
                family = "binomial", data = qpcr, na.action = na.omit)) %>% 
  map_df(.f = ~tidy(.x,  conf.int = TRUE)) %>%
  filter(term != "(Intercept)")

death_stats %>% filter(p.value < 0.05)
# A tibble: 0 × 7
# … with 7 variables: term <chr>, estimate <dbl>, std.error <dbl>, statistic <dbl>,
#   p.value <dbl>, conf.low <dbl>, conf.high <dbl>

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
multi_resistant_stats <-  glm(resistance ~., data = data_resistant,
                              family="binomial") %>% tidy(conf.int = TRUE)
multi_resistant_stats %>% filter(p.value < 0.05)
# A tibble: 0 × 7

data_death <- qpcr %>% select(-recurrence, -resistance) 
multi_death_stats <-  glm(death ~., data = data_death,
                              family="binomial") %>% tidy(conf.int = TRUE)
multi_death_stats %>% filter(p.value < 0.05)
# A tibble: 0 × 7

data_recurrence <- qpcr %>% select(-death, -resistance) 
multi_recurrence_stats <-  glm(recurrence ~., data = data_recurrence,
                          family="binomial") %>% tidy(conf.int = TRUE)
multi_recurrence_stats %>% filter(p.value < 0.05)
# A tibble: 0 × 7
