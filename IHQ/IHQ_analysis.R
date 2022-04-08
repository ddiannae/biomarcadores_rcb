library(openxlsx)
library(stringi)
library(stringr)
library(limma)
library(affy) ## falta
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38) ## falta
library(org.Hs.eg.db)
library(clusterProfiler)
library(Biobase)
library(AnnoProbe)
library(jtools)
library(ranger)
library(caret)
library(varSelRF)
library(Hmisc)
library(ggpubr)
library(gridExtra)
library(corrplot)
library(MKmisc)
library(pROC)
library(plotly)
library(tidymodels)
library(survival)
library(lubridate)
library(survminer)
library(patchwork)


IHQ.df <- read.csv("data_IHQ/IHQ_2022.csv", header = TRUE)
##### Modelos de regression logistica para predecir respuesta con IHQs
##### Preparamos el data frame los distintos analisis



covariates <- c("Age","BMI", "cT","cN","TNM","Tum_size", 
                "Histology","Menopause",
                "Ki67","Her2")


IHQ.df$Response <- factor(ifelse(IHQ.df$Resistant == 1, "Resistant", "Sensitive"), levels = c("Sensitive","Resistant"))
IHQ.df$Menopause <- factor(ifelse(IHQ.df$Menop_Status == 1, "Positive", "Negative"),
                           levels = c("Negative","Positive"))
IHQ.df$Her2 <- factor(ifelse(IHQ.df$HER2 == 1, "Positive", "Negative"),
                      levels = c("Negative","Positive"))

### Aqui separamos los que tienen GRIApreNAC en score Low y High en base
### a la media. Si crees que es mejor usar otro punto de corte, 
### es mejor cambiarlo porque eso puede afectar los analsis de supervivencia
IHQ.df$IHC_class_GRIApreNAC <- factor(ifelse(IHQ.df$IHC_score_GRIApreNAC <= median(IHQ.df$IHC_score_GRIApreNAC, na.rm=T), "Low", "High"), 
                                      levels = c("Low","High"))
IHQ.df$IHC_class_NDUFpreNAC <- factor(ifelse(IHQ.df$IHC_score_NDUFpreNAC <= median(IHQ.df$IHC_score_NDUFpreNAC, na.rm=T), "Low", "High"), 
                                      levels = c("Low","High"))


IHQ.df$Last_fup_Date <- ymd(IHQ.df$Last_fllowup)
IHQ.df$Admission_Date <- ymd(IHQ.df$Admission)
IHQ.df$Treatment_Start_Date <- ymd(IHQ.df$NAC_start)
IHQ.df$Recurrance_Date <- ymd(IHQ.df$Recurrence_date)
IHQ.df$Death_Date <- ymd(IHQ.df$Death_date)
IHQ.df$Death.Status <- ifelse(is.na(IHQ.df$Death_Date), 1, 2)
IHQ.df$Recurrance.Status <- ifelse(is.na(IHQ.df$Recurrance_Date), 1, 2)

IHQ.df$os.yrs <- as.duration(IHQ.df$Treatment_Start_Date %--% IHQ.df$Last_fup_Date) / dyears(1)
IHQ.df$censor.time <- as.duration(IHQ.df$Treatment_Start_Date %--% IHQ.df$Last_fup_Date) / dyears(1)
IHQ.df$recurrance.time <- as.duration(IHQ.df$Treatment_Start_Date %--% IHQ.df$Recurrance_Date) / dyears(1)
IHQ.df$dfs.yrs <- pmin(IHQ.df$censor.time,IHQ.df$recurrance.time, na.rm = TRUE)

### Analisis a y b - Estadistica simple y boxplots para
### ver diferencias entre a) preNAC y posNAC en NDUF
### y entre sensibles y resistentes en preNAC de NDUF y GRIA
### Filtre los que tenian IHC_score de 0, pero no se si es correcto hacerlo

IHQ.res.df <- filter(IHQ.df, !is.na(Response), PR == 1, !is.na(IHC_score_NDUFpreNAC))


nduf.pre <- filter(IHQ.df,!is.na(IHC_score_NDUFpreNAC), Response == "Resistant")$IHC_score_NDUFpreNAC
nduf.pos <- filter(IHQ.df,!is.na(IHC_score_NDUFposNAC), Response == "Resistant")$IHC_score_NDUFposNAC

IHQ.nduf.df <- data.frame(IHS_score = c(nduf.pre,nduf.pos),
                          sample = c(rep("PreNAC",length(nduf.pre)),
                                     rep("PosNAC",length(nduf.pos))),
                          stringsAsFactors = FALSE)

IHQ.nduf.df$sample <- factor(IHQ.nduf.df$sample, levels = c("PreNAC","PosNAC"))

filter(IHQ.nduf.df, IHS_score > 0) %>%
  ggplot(aes(x = sample, y = IHS_score)) +
  geom_boxplot(width = 0.5, lwd = 1.5, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.1, size = 2, aes(color = sample)) +
  xlab("") +
  ylab("NDUF IHS score") +
  ylim(0,15) +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test")


nduf.pre <- filter(IHQ.df,!is.na(IHC_score_NDUFpreNAC))[,c("IHC_score_NDUFpreNAC","Response")]
colnames(nduf.pre) <- c("IHC_score_preNAC","Response")
nduf.pre$gene <- "NDUF"
gria.pre <- filter(IHQ.df,!is.na(IHC_score_GRIApreNAC))[,c("IHC_score_GRIApreNAC","Response")]
colnames(gria.pre) <- c("IHC_score_preNAC","Response")
gria.pre$gene <- "GRIA"
IHQ.pre.df <- rbind(nduf.pre, gria.pre)


filter(IHQ.pre.df, IHC_score_preNAC > 0) %>%
  ggplot(aes(x = Response, y = IHC_score_preNAC)) +
  geom_boxplot(width = 0.5, lwd = 1.5, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.1, size = 2, aes(color = Response)) +
  xlab("") +
  ylab("IHS score preNAC") +
  ylim(0,15) +
  scale_color_manual(values = c("darkblue","orange")) +
  theme_bw(base_size = 18) +
  stat_compare_means(method = "wilcox.test",label.y = 13) +
  facet_wrap(~gene, nrow = 1)

### Analisis 1 - Queremos ver que variables se correlacionan mas con
### los scores de IHQ de nuestros genes, para no incluirlos en los modelos
### de regresion, ya que esto puede ocacionar problemas

## Opcion 1: Regresiones logisticas univariables -
## Esto es mejor cuando al menos una de las dos variables a comparar es continua
## Tambien puede ser buena historia para reportar si hay una correlacion alta

## Las variables dependientes tiene que ser binarias, asi que en este caso
## dicotomizamos el IHC score and "Low" y "High" (tu sabras cual es el mejor punto de corte
## que yo use solo la media)

dependant.variables <- c("IHC_class_GRIApreNAC","IHC_class_NDUFpreNAC")
other.variables <- c("Age","BMI", "cT","cN","TNM","Tum_size", 
                     "Histology","Menopause",
                     "Ki67","Her2","Response")

IHQ.logreg.uni.assoc <- list()
IHQ.logreg.uni.assoc.df <- c()
for (dv in dependant.variables) {
  for (ov in other.variables) {
    tmp.df <- IHQ.res.df[!is.na(IHQ.res.df[,dv]) & !is.na(IHQ.res.df[,ov]),]
    dp.f <- paste(dv, "~", ov)
    glm.fit <- glm(as.formula(dp.f), data = tmp.df, 
                   family = "binomial")
    IHQ.logreg.uni.assoc[[dp.f]] <-  glm.fit
    t <- summ(glm.fit)
    tmp.df <- data.frame(dep.var = dv, other.var = ov, 
                         Estimate = t$coeftable[2,1],
                         Standard.Error = t$coeftable[2,2],
                         Zvalue = t$coeftable[2,3],
                         pval = t$coeftable[2,4])
    IHQ.logreg.uni.assoc.df <- rbind(IHQ.logreg.uni.assoc.df, tmp.df)
  }
}
IHQ.logreg.uni.assoc.df$padj <- p.adjust(IHQ.logreg.uni.assoc.df$pval, method = "BH")

IHQ.logreg.uni.assoc.df <- arrange(IHQ.logreg.uni.assoc.df, pval, -Estimate) %>% as.data.frame()

# Podemos aqui checar que variables salieron 
# significantes. En este caso no saluo ninguna
# asi que pongo el punto de corte muy bajo nada mas para que veas
filter(IHQ.logreg.uni.assoc.df, pval < 0.1)

# Y aqui puedes ver el resumen del modelo univariable
# para predecir IHC_class_GRIApreNAC con Tum_siz1
summary(IHQ.logreg.uni.assoc$`IHC_class_GRIApreNAC ~ Tum_size`)

## Opcion 2: Fisher test -
## Esto es mejor cuando al menos una de las dos variables a comparar son
## discretas (a lo mucho 4 grupos en cada variable)

## Aqui solo pongo de ejemplo GRIA vs Menopousia que no son
## significantes, si quieres hacer mas tienes que dicotomizar el resto

discrete.vars <- c("IHC_class_GRIApreNAC","IHC_class_NDUFpreNAC",
                   "Menopause","Histology","Response")

IHQ.fisher.test.df <- c()
for (i in 1:(length(discrete.vars) - 1)) {
  vr1 <- discrete.vars[i]
  for (j in (i+1):length(discrete.vars)) {
    vr2 <- discrete.vars[j]
    ft <- fisher.test(table(IHQ.res.df[,vr1], IHQ.res.df[,vr2]))
    tmp.df <- data.frame(var1 = vr1, var2 = vr2,
                         pval = ft$p.value,
                         odds = ft$estimate,
                         logoddss = log2(ft$estimate),
                         lower.conf.int = ft$conf.int[1],
                         upper.conf.int = ft$conf.int[2])
    IHQ.fisher.test.df <- rbind(IHQ.fisher.test.df, tmp.df)
  }
}
IHQ.fisher.test.df <- arrange(IHQ.fisher.test.df, pval, -odds) %>% as.data.frame()

IHQ.fisher.test.df$padj <- p.adjust(IHQ.fisher.test.df$pval, method = "BH")


grid.table(IHQ.fisher.test.df)

## Opcion 3: Correlaciones -
## Esto es mejor cuando al menos una de las dos variables a comparar son
## discretas (a lo mucho 4 grupos en casa variable)

IHQ.cor <- rcorr(as.matrix(IHQ.res.df[,c("IHC_score_GRIApreNAC","IHC_score_NDUFpreNAC","Age",
                                         "BMI","Ki67","cT","cN","TNM","Tum_size",
                                         "RCB_class")]))

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

IHQ.cor.df <- flattenCorrMatrix(IHQ.cor$r, IHQ.cor$P)
IHQ.cor.df <- arrange(IHQ.cor.df, p)
IHQ.cor.df$padj <- p.adjust(IHQ.cor.df$p, method = "BH")

filter(IHQ.cor.df, p < 0.05)
filter(IHQ.cor.df, padj < 0.05)

IHQ.cor.plt <- cor(as.matrix(IHQ.res.df[,c("IHC_score_GRIApreNAC","IHC_score_NDUFpreNAC","Age",
                                           "BMI","Ki67","cT","cN","TNM","Tum_size",
                                           "RCB_class")]), use = "pairwise.complete.obs")

corrplot(IHQ.cor.plt, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)


### Conclusiones de este analisis: Parece GRIA y Ki67 estan correlacionados
### Por lo que tal vez sea mejor quitar Ki67 del modelos logisticos lineales
### al evaluar GRIA para ver cual de los dos es mejor predictor de respuesta
### De igual forma cT, Tum_size y TNM parecen algo correlacionados
### asi que igual mejor sea solo escoger una de esas para incluir en el modelo

  
  
#### ANALISIS 2 - Ver como se relacionan los scores de IHQ con respuesta
#### con y sin ajuste a covariables
  
### Primero hacemos modelos univariables

dependant.variable <- "Response"
independant.variables <- c("IHC_score_GRIApreNAC","IHC_score_NDUFpreNAC")
covariates <- c("Age","BMI", "cT","cN","TNM","Tum_size", 
                "Histology","Menopause",
                "Ki67","Her2")


IHQ.response.univariate.list <- list()
IHQ.response.univariate.df <- c()
for (iv in c(covariates, independant.variables)) {
  tmp.df <- IHQ.res.df[!is.na(IHQ.res.df[,iv]),]
  uni.f <- paste(dependant.variable, "~", iv)
  glm.fit <- IHQ.response.univariate.list[[iv]] <- glm(as.formula(uni.f), data = tmp.df, 
                                                       family = "binomial")
  IHQ.response.univariate.list[[iv]] <- glm.fit
  t <- summ(glm.fit)
  tmp.df <- data.frame(dep.var = "Response", other.var = iv, 
                       Estimate = t$coeftable[2,1],
                       Standard.Error = t$coeftable[2,2],
                       Zvalue = t$coeftable[2,3],
                       pval = t$coeftable[2,4])
  IHQ.response.univariate.df <- rbind(IHQ.response.univariate.df, tmp.df)
}
IHQ.response.univariate.df$padj <- p.adjust(IHQ.response.univariate.df$pval, method = "BH")
IHQ.response.univariate.df <- arrange(IHQ.response.univariate.df,pval)

filter(IHQ.response.univariate.df, pval < 0.1)

summary(IHQ.response.univariate.list$IHC_score_GRIApreNAC)


# Ahora hacemos Modelos multiiveriables

dependant.variable <- "Response"
independant.variables <- c("IHC_score_GRIApreNAC","IHC_score_NDUFpreNAC")

# Aqui pongo de ejempo estas covariables para incluir
# en el modelo completo, pero puedes poner los que creas
# necesarios.  
covariates <- c("Age","BMI","Menopause") 

IHQ.response.multivariate.list <- list()
for (iv in independant.variables) {
  dp.f <- paste(c(covariates,iv), collapse = " + ")
  multi.f <- paste(dependant.variable, "~", dp.f)
  IHQ.response.multivariate.list[[iv]] <- glm(as.formula(multi.f), data = IHQ.res.df, 
                                               family = "binomial")
}

summary(IHQ.response.multivariate.list$IHC_score_NDUFpreNAC)



############ DIANA HASTA AQUI
#### Analisis 4 Random forest para ver cuales
#### variables son mas importantes en el modelo completo (usando todas las variables)
#### en predecir respuesta.
#### El objetivo no es tanto hacer una buena prediccion sino evaluar
#### cuales variables son mas robustas y son mas importantes


dependant.variable <- "Response"
all.features <- c("Age","BMI","Tum_size","cT","cN","TNM","Ki67","Histology",
                  "IHC_score_GRIApreNAC","IHC_score_NDUFpreNAC","Her2")

IHQ.res.df$ID <- paste0("IHQ_",1:nrow(IHQ.res.df))
rownames(IHQ.res.df) <- IHQ.res.df$ID  
ihq.res <- filter(IHQ.res.df, Response == "Resistant")$ID
ihq.sns <- filter(IHQ.res.df, Response == "Sensitive")$ID

### Aqui preparamos las variables para random forest, deben
### estar todas centradas y escaladas

tmp.rf.df <- IHQ.res.df[,c("Response",all.features)]
tmp.rf.df$Age <- as.numeric(scale(tmp.rf.df$Age, center = TRUE, scale = TRUE))
tmp.rf.df$BMI <- as.numeric(scale(tmp.rf.df$BMI, center = TRUE, scale = TRUE))
tmp.rf.df$Tum_size <- as.numeric(scale(tmp.rf.df$Tum_size, center = TRUE, scale = TRUE))
tmp.rf.df$cT <- as.numeric(scale(tmp.rf.df$cT, center = TRUE, scale = TRUE))
tmp.rf.df$cN <- as.numeric(scale(tmp.rf.df$cN, center = TRUE, scale = TRUE))
tmp.rf.df$TNM <- as.numeric(scale(tmp.rf.df$TNM, center = TRUE, scale = TRUE))
tmp.rf.df$Ki67 <- as.numeric(scale(tmp.rf.df$Ki67, center = TRUE, scale = TRUE))
tmp.rf.df$Histology <- as.numeric(scale(tmp.rf.df$Histology, center = TRUE, scale = TRUE))
tmp.rf.df$IHC_score_GRIApreNAC <- as.numeric(scale(tmp.rf.df$IHC_score_GRIApreNAC, center = TRUE, scale = TRUE))
tmp.rf.df$IHC_score_NDUFpreNAC <- as.numeric(scale(tmp.rf.df$IHC_score_NDUFpreNAC, center = TRUE, scale = TRUE))

dp.f <- paste0(all.features, collapse = " + ")

#### Para que este balanceadas las clases
all.res.comb <- sample(combn(ihq.res, 8, simplify = FALSE), size = 200, replace = FALSE)

fold <- 1
IHQ.var.importances.df <- c()
for (sns in ihq.sns) {
  for (i in 1:length(all.res.comb)) {
    train.smps <- c(setdiff(ihq.sns,sns),all.res.comb[[i]])
    rf.rang.obj <- ranger(formula = as.formula(paste(dependant.variable, "~",  dp.f)), 
                          data = tmp.rf.df[train.smps,], num.trees = 10000, num.threads = 2,
                          importance = "impurity")
    rf.imp <- rf.rang.obj$variable.importance
    tmp.df <- data.frame(variable = names(rf.imp),
                         importance = rf.imp,
                         fold = fold, stringsAsFactors = TRUE)
    IHQ.var.importances.df <- rbind(IHQ.var.importances.df, tmp.df)
    fold <- fold + 1
  }
}

IHQ.var.mean.df <- IHQ.var.importances.df %>% group_by(variable) %>%
  summarise(mean.importance = mean(importance),
            mean.variance = var(importance)) %>% arrange(-mean.importance) %>%
  as.data.frame()
IHQ.var.importances.df$variable <- factor(IHQ.var.importances.df$variable,
                                          levels = IHQ.var.mean.df$variable)

### Graficamos las importancias y
### podemos ver las mas importantes
ggplot(IHQ.var.importances.df, 
       aes(x = variable, y = importance, fill = variable)) +
  geom_boxplot(width = 0.4) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

#### Analisis 5 Modelos de regresion logistica vs Random forest 
#### univariables para ver cuales predicen mejor Respuesta.
#### Hace validacion cruzada de 3-folds y repetida 10 veces
#### hacer estadistica y ver cuales variables son mas robustas
#### para predecir.

folds <- 3
all.features <- c("Age","BMI","Tum_size","cT","cN","TNM","Ki67",
                  "IHC_score_GRIApreNAC","IHC_score_NDUFpreNAC")

cvIndex <- createMultiFolds(IHQ.res.df$Response, 
                            k = folds, times = 10)

IHQ.cv.scores.df <- c()
IHQ.cv.prcurves.df <- c()
IHQ.cv.aucurves.df <- c()

for (iv in all.features) {
  tmp.rf2.df <- tmp.rf.df[!is.na(tmp.rf.df[,iv]),]
  tmp.lr.df <- IHQ.res.df[!is.na(IHQ.res.df[,iv]),]
  for (fd in names(cvIndex)) {
    train.smps <- cvIndex[[fd]]
    dp.f <- paste(dependant.variable, "~", iv)
    
    glm.fit <- glm(as.formula(dp.f), data = tmp.lr.df[train.smps,], 
                   family = "binomial")
    
    glm.pred <- predict(glm.fit,tmp.lr.df[-train.smps,iv,drop = FALSE], 
                        type = "response")
    
    
    preds.df <- data.frame(truth = tmp.lr.df[-train.smps,dependant.variable],
                           Sensitive = 1 - glm.pred,
                           Resistant = glm.pred,
                           predicted = factor(ifelse(glm.pred > 0.5, "Resistant","Sensitive"),
                                              levels = c("Resistant","Sensitive")))
    
    perf.df <- perfMeasures(pred = preds.df$Resistant, truth = tmp.lr.df[-train.smps,dependant.variable], namePos = "Resistant")
    perf.df <- filter(perf.df, Measure %in% c("F1 score","accuracy (ACC)"))
    
    
    glm.roc <- roc_auc(preds.df, truth, Resistant)$.estimate
    glm.roc <- max(glm.roc, 1 - glm.roc)
    glm.prc <- pr_auc(preds.df, truth, Resistant)$.estimate
    glm.prc <- max(glm.prc, 1 - glm.prc)
    
    glm.tmp.df <- data.frame(IV = iv, fold = fd,
                             Measure = c(perf.df$Measure, "auROC","auPRC"),
                             Value = c(perf.df$Value, 
                                       glm.roc,
                                       glm.prc),
                             model = "LR",
                             stringsAsFactors = FALSE)
    
    IHQ.cv.scores.df <- rbind(IHQ.cv.scores.df, glm.tmp.df)
    
    prc.curve.df <- as.data.frame(pr_curve(preds.df, truth, Resistant))
    
    colnames(prc.curve.df) <- c("cutoff","Recall", "Precision")
    prc.curve.df$IV <- iv
    prc.curve.df$fold <- fd
    prc.curve.df$model <- "LR"
    IHQ.cv.prcurves.df <- rbind(IHQ.cv.prcurves.df, prc.curve.df)
    
    auc.curve.df <- as.data.frame(roc_curve(preds.df, truth, Resistant))
    
    colnames(auc.curve.df) <- c("cutoff","Recall", "Precision")
    auc.curve.df$IV <- iv
    auc.curve.df$fold <- fd
    auc.curve.df$model <- "LR"
    IHQ.cv.aucurves.df <- rbind(IHQ.cv.aucurves.df, auc.curve.df)
    
    rf.rang.obj <- ranger(formula = as.formula(dp.f), 
                          data = tmp.rf2.df[train.smps,], 
                          num.trees = 10000, num.threads = 2,
                          importance = "impurity", probability = TRUE)
    
    preds <- predict(rf.rang.obj,tmp.rf2.df[-train.smps,iv,drop = FALSE], 
                     type = "response")
    
    
    preds.df <- data.frame(truth = tmp.rf2.df[-train.smps,dependant.variable],
                           Sensitive = preds$predictions[,"Sensitive"],
                           Resistant = preds$predictions[,"Resistant"],
                           predicted = factor(ifelse(preds$predictions[,"Resistant"] > 0.5, "Resistant","Sensitive"),
                                              levels = c("Resistant","Sensitive")))
    
    perf.df <- perfMeasures(pred = preds$predictions[,2], truth = tmp.rf2.df[-train.smps,dependant.variable], namePos = "Resistant")
    perf.df <- filter(perf.df, Measure %in% c("F1 score","accuracy (ACC)"))
    
    
    rf.roc <- roc_auc(preds.df, truth, Resistant)$.estimate
    rf.roc <- max(rf.roc, 1 - rf.roc)
    rf.prc <- pr_auc(preds.df, truth, Resistant)$.estimate
    rf.prc <- max(rf.prc, 1 - rf.prc)
    
    
    tmp.df <- data.frame(IV = iv, fold = fd,
                         Measure = c(perf.df$Measure, "auROC","auPRC"),
                         Value = c(perf.df$Value, rf.roc,
                                   rf.prc),
                         model = "RF",
                         stringsAsFactors = FALSE)
    
    IHQ.cv.scores.df <- rbind(IHQ.cv.scores.df, tmp.df)
    
    prc.curve.df <- as.data.frame(pr_curve(preds.df, truth, Resistant))
    
    colnames(prc.curve.df) <- c("cutoff","Recall", "Precision")
    prc.curve.df$IV <- iv
    prc.curve.df$fold <- fd
    prc.curve.df$model <- "RF"
    IHQ.cv.prcurves.df <- rbind(IHQ.cv.prcurves.df, prc.curve.df)
    
    auc.curve.df <- as.data.frame(roc_curve(preds.df, truth, Resistant))
    
    colnames(auc.curve.df) <- c("cutoff","Recall", "Precision")
    auc.curve.df$IV <- iv
    auc.curve.df$fold <- fd
    auc.curve.df$model <- "RF"
    IHQ.cv.aucurves.df <- rbind(IHQ.cv.aucurves.df, auc.curve.df)
  }
}

IHQ.cv.scores.df$IV <- factor(IHQ.cv.scores.df$IV,levels = all.features)

# Aqui graficamos los scores de prediccion para los modelos
# logisticos (LR) y los de random forest (RF)
# Como las clases esta desbalanaceadas, es mejor
# hacerle caso al area under the precios-recall cruve (auPRC)
# o al F1-score
filter(IHQ.cv.scores.df) %>%
  ggplot(aes(x = IV, y = Value, fill = model)) +
  geom_boxplot(width = 0.4) +
  ylab("Performance") +
  theme_bw(base_size = 18) +
  xlab("Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Measure, nrow = 2)



#### Analisis 6 Aqui hacemos supervivencia para
#### overall survival (OS) y disease-free survival (DFS)

# Analisis univariable: Primero hacemos curvas kaplan-meyer para ver
# si pacientes con nivel alto/bajo GRIA o NDUF si tienen significativamente
# mas o menor probabilidad se sobrevivir o recurrir despues de la terapia
os.GRIApreNAC.fit <- survfit(Surv(os.yrs, Death.Status) ~ IHC_class_GRIApreNAC, 
                             data = filter(IHQ.df, !is.na(IHC_class_GRIApreNAC)))

ggsurvplot(os.GRIApreNAC.fit, data = filter(IHQ.df, !is.na(IHC_class_GRIApreNAC)),
           size = 1, palette = c("Low" = "blue", "High" = "red"), 
           legend.labs = c("Low","High"),
           pval = TRUE,  
           ggtheme = theme_bw(base_size = 18)) + guides(color = guide_legend(title = "GRIA pre-NAC"))

os.NDUFpreNAC.fit <- survfit(Surv(os.yrs, Death.Status) ~ IHC_class_NDUFpreNAC, 
                             data = filter(IHQ.df, !is.na(IHC_class_NDUFpreNAC)))

ggsurvplot(os.NDUFpreNAC.fit, data = filter(IHQ.df, !is.na(IHC_class_NDUFpreNAC)),
           size = 1, palette = c("Low" = "blue", "High" = "red"), 
           legend.labs = c("Low","High"),
           pval = TRUE,  
           ggtheme = theme_bw(base_size = 18)) + guides(color = guide_legend(title = "NDUF pre-NAC"))


dfs.GRIApreNAC.fit <- survfit(Surv(dfs.yrs, Death.Status) ~ IHC_class_GRIApreNAC, 
                              data = filter(IHQ.df, !is.na(IHC_class_GRIApreNAC)))

ggsurvplot(dfs.GRIApreNAC.fit, data = filter(IHQ.df, !is.na(IHC_class_GRIApreNAC)),
           size = 1, palette = c("Low" = "blue", "High" = "red"), 
           legend.labs = c("Low","High"),
           pval = TRUE,  
           ggtheme = theme_bw(base_size = 18)) + 
  guides(color = guide_legend(title = "GRIA pre-NAC")) + ylab("Recurrance probability")

dfs.NDUFpreNAC.fit <- survfit(Surv(dfs.yrs, Death.Status) ~ IHC_class_NDUFpreNAC, 
                              data = filter(IHQ.df, !is.na(IHC_class_NDUFpreNAC)))

ggsurvplot(dfs.NDUFpreNAC.fit, data = filter(IHQ.df, !is.na(IHC_class_NDUFpreNAC)),
           size = 1, palette = c("Low" = "blue", "High" = "red"), 
           legend.labs = c("Low","High"),
           pval = TRUE,  
           ggtheme = theme_bw(base_size = 18)) + 
  guides(color = guide_legend(title = "NDUF pre-NAC")) +
  ylab("Recurrance probability")



### Podemos hacer lo mismo para todas las variables
### y ver cuales son mas correlacionadas
### con supervivencia usando cox-hazard
independant.variables <- c("IHC_score_GRIApreNAC","IHC_score_NDUFpreNAC")
covariates <- c("Age","BMI", "cT","cN","TNM","Tum_size", 
                "Histology","Menopause",
                "Ki67","Her2")


IHQ.coxhazard.univariate.os.list <- list()
IHQ.coxhazard.univariate.dfs.list <- list()
IHQ.coxhazard.univariate.df <- c()
for (iv in c(covariates, independant.variables)) {
  IHQ.tmp.df <- IHQ.df[!is.na(IHQ.df[,iv]),]
  
  cox.f <- as.formula(paste("Surv(os.yrs, Death.Status)", "~", iv))
  cox.fit <-  coxph(cox.f, data = IHQ.tmp.df)
  IHQ.coxhazard.univariate.os.list[[iv]] <- cox.fit
  t <- summary(cox.fit)
  
  tmp.df <- data.frame(Analysis = "OS", other.var = iv, 
                       Estimate = t$coefficients[1],
                       Standard.Error = t$coefficients[3],
                       Zvalue = t$coefficients[4],
                       pval = t$coefficients[5])
  IHQ.coxhazard.univariate.df <- rbind(IHQ.coxhazard.univariate.df, tmp.df)
  
  cox.f <- as.formula(paste("Surv(dfs.yrs, Recurrance.Status)", "~", iv))
  cox.fit <-  coxph(cox.f, data = IHQ.tmp.df)
  IHQ.coxhazard.univariate.dfs.list[[iv]] <- cox.fit
  t <- summary(cox.fit)
  
  tmp.df <- data.frame(Analysis = "DFS", other.var = iv, 
                       Estimate = t$coefficients[1],
                       Standard.Error = t$coefficients[3],
                       Zvalue = t$coefficients[4],
                       pval = t$coefficients[5])
  IHQ.coxhazard.univariate.df <- rbind(IHQ.coxhazard.univariate.df, tmp.df)
}

IHQ.coxhazard.univariate.df$Analysis <- factor(IHQ.coxhazard.univariate.df$Analysis,
                                               levels = c("OS","DFS"))

IHQ.coxhazard.univariate.df$padj <- p.adjust(IHQ.coxhazard.univariate.df$pval, method = "BH")
IHQ.coxhazard.univariate.df <- IHQ.coxhazard.univariate.df %>% 
  arrange(Analysis, pval)

filter(IHQ.coxhazard.univariate.df, Analysis == "OS", pval < 0.1)

filter(IHQ.coxhazard.univariate.df, Analysis == "DFS", pval < 0.1)

summary(IHQ.coxhazard.univariate.os.list$IHC_score_GRIApreNAC)


### Ahora hacemos Modelos multiiveriables con 
### cox regression para ver riesgos


independant.variables <- c("IHC_score_GRIApreNAC","IHC_score_NDUFpreNAC")

covariates1 <- c("Age","Her2","Ki67") # aqui agrega las covariables que creas necesarias
covariates2 <- c("Age","Her2","Menopause","Tumor_size")


IHQ.response.multivariate.os.list <- list()
for (iv in independant.variables) {
  multi.f <- as.formula(paste("Surv(os.yrs, Death.Status)", "~", paste(c(covariates1,iv), collapse = "+")))
  IHQ.response.multivariate.os.list[[iv]] <- coxph(multi.f, data = IHQ.df)
}

IHQ.response.multivariate.dfs.list <- list()
for (iv in independant.variables) {
  multi.f <- as.formula(paste("Surv(dfs.yrs, Recurrance.Status)", "~", paste(c(covariates2,iv), collapse = "+")))
  IHQ.response.multivariate.dfs.list[[iv]] <- coxph(multi.f, data = IHQ.df)
}


summary(IHQ.response.multivariate.os.list$IHC_score_NDUFpreNAC)
