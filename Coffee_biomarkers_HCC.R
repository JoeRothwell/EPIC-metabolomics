# Coffee biomarker data from Profinder for HCC study. Read in data generated from Profinder
# All samples, blanks and QCs

# Data preparation ----
library(tidyverse)
coffee4hcc <- function(){
  cycloproval <- read_csv("data/Cyclo pro val HCC data_1.csv") %>% slice(1)

  #pos mode data
  pos <- read_csv("Coffee 8 pos biomarkers HCC data.csv") %>% bind_rows(cycloproval) %>% 
    select("Compound Name", contains("Area")) %>% t
  colnames(pos) <- pos[1, ]
  pos           <- pos[-1, ]
  pos.samp      <- pos[1:258, ]

  #neg mode data
  neg <- read_csv("Coffee 2 neg biomarkers HCC data.csv") %>% select("Compound Name", contains("Area")) %>% t
  colnames(neg) <- neg[1, ]
  neg           <- neg[-1, ]
  neg.samp      <- neg[1:258, ]

  #get labID (same for pos and neg), first 258 observations are controls
  sampnames <- rownames(neg)[1:258]
  LabIDs    <- data.frame(ID = sampnames) %>% 
  separate(ID, into = c("other", "other2", "LabID", "No"), sep=c(31,41,51))

  #add LabID as rownames and merge pos and neg biomarkers
  all.samp           <- cbind(LabID = LabIDs$LabID, pos.samp, neg.samp)
  rownames(all.samp) <- NULL
  all.samp           <- as.tibble(all.samp)
  #all.samp2          <- lapply(all.samp, FUN=as.numeric) %>% data.frame

  #Get metadata. Read in from uncompressed SAS file and subset subject metadata
  library(haven)
  hcc <- read_sas("data/merged_untarg1.sas7bdat") %>% 
    select(Id_Bma, Idepic_Samp:alc_drinker_m) %>%
    separate("Id_Bma", into = c("Date", "LabID", "Mode"), sep=c(9,19))
  
  #For liver function
  lf <- read_sas("data/live_caco1.sas7bdat") %>% 
    filter(Match_Caseset != 387 | Idepic != "22____22304249") %>%
    select(Idepic, Liver_Function_Score, Liver_Function_Score_C)

  #For coffee intakes
  cof <- readRDS("prepdata/Coffee intakes all EPIC.rds")
  
  hcc <- hcc %>% left_join(cof, by="Idepic")
  hcc <- hcc %>% left_join(lf, by="Idepic")

  #join coffee biomarker data
  chcc <- inner_join(hcc, all.samp, by="LabID") %>% select(-Date, -Mode) %>%
    #rename compound names for easier analysis in STATA and convert to numeric
  rename(Trigonelline       = `Trigonelline (3-Carboxy-1-methylpyridinium betaine)`,
         Quinic_acid        = `Quinic acid`, 
         Hippuric_acid      = `Hippuric acid`, 
         Cyclo_leu_pro      = `Cyclo(leucyl-prolyl)`,
         Cyclo_isol_pro     = `Cyclo(isoleucyl-prolyl)`,
         Cyclo_pro_val      = `Cyclo(prolyl-valyl)`,
         Cat_sulf           = `Catechol sulfate`,
         AAMU               = `5-Acetylamino-6-amino-3-methyluracil (AAMU)`) %>% 
         mutate_at(vars(Trigonelline:AAMU), as.numeric)

  #Make intensity matrix
  ints <- chcc %>% select(Trigonelline:AAMU) %>% as.matrix

  #Convert numeric columns to factor
  chcc$Smoke_Stat    <- as.factor(chcc$Smoke_Stat)
  chcc$alc_drinker_m <- as.factor(chcc$alc_drinker_m)
  chcc$L_School      <- as.factor(chcc$L_School)
  #chcc$Caselive_Crs  <- as.factor(chcc$Caselive_Crs)
  
  return(chcc)
  #saveRDS(coffee.hcc, "HCC coffee biomarkers and metadata.rds")
  #write to .dta (Stata)
  #library(foreign)
  #write.dta(coffee.hcc, "D://HCC_coffee_metabolites1.dta")
}
chcc <- coffee4hcc()

# saveRDS(chcc, "HCC coffee biomarkers and metadata.rds")
chcc <- readRDS("prepdata/HCC coffee biomarkers and metadata.rds")


# Association HCC coffee ----


# Boxplot cases-controls
ggplot(chcc, aes(x=as.factor(Caselive_Crs), y=sqrt(QGE130301))) + 
  geom_boxplot(fill = "grey", outlier.colour = "white") + 
  theme_classic() + geom_jitter(width=0.1, shape=1)

# Matched boxplot
chcc <- chcc %>% mutate(Case = ifelse(Caselive_Crs == 0, -3, 3), 
                        Case2 = ifelse(Caselive_Crs == 0, -1.5, 1.5))
ggplot(chcc, aes(x = Case, y= sqrt(QGE130301), group = Caselive_Crs)) +
  geom_boxplot(fill="dodgerblue", width=2) + theme_bw() +
  geom_line(aes(x = Case2, group=Match_Caseset), colour="grey") +
  xlab("HCC controls vs cases") + ylab("sqrt(coffee intake)") +
  geom_point(aes(x = Case2))

# Replicate categories of Aleksandrova's paper
chcc$cof.cat <- cut(chcc$QGE130301, breaks = c(0, 300, 450, 600, 1600), labels=c(1:4), include.lowest = T)
table(chcc$cof.cat, chcc$Caselive_Crs)

# Coloured effects plot
ggplot(chcc, aes(x = as.factor(Caselive_Crs), y= sqrt(QGE130301))) +
  geom_line(aes(group=Match_Caseset, colour=as.factor(cof.cat))) + theme_classic() +
  geom_point() +
  xlab("Controls vs cases")

# Linear model
fit <- lm(log(chcc$QGE130301 + 1) ~ chcc$Caselive_Crs)

library(survival)
cof.CLR <- clogit(Caselive_Crs ~ log2(QGE130301+1) + Bmi_C + Alc_Re + Smoke_Stat + 
        alc_drinker_m + L_School + Pa_Mets + Waist_C + strata(Match_Caseset), data = chcc)


# Association HCC coffee metabolites ----
# check missings and equal cases and controls
library(Amelia)
missmap(chcc, rank.order=F)
table(chcc$Caselive_Crs)

# Conditional logistic model

# See email from Magda "metabo" 16-feb-17
# Covariates are Bmi_C, Alc_Re, Smoke_Stat, alc_drinker_m, L_School, Pa_Mets, Waist_C
# Function to fit a CLR for each metabolite

CLR.coffee <- function(adj = T) {
  
  #Models for 10 biomarkers, raw or adjusted model
  library(survival)
  fit.CLR <- if(adj == T) { 
    function(x) {
            clogit(Caselive_Crs ~ log2(x) + #log2(QGE130301+1) + 
            Bmi_C + Alc_Re + Smoke_Stat + alc_drinker_m + 
            L_School + Pa_Mets + Waist_C + strata(Match_Caseset), data = chcc)
            }
                } else { 
    function(x) {  
            clogit(Caselive_Crs ~ log2(x) + strata(Match_Caseset), data = chcc)
            }
      }

  #Make intensity matrix
  library(tidyverse)
  ints <- chcc %>% select(Trigonelline:AAMU) %>% as.matrix
  #adjusted or unadjusted models
  mods <- apply(ints, 2, fit.CLR)

  Pcor  <- unlist(sapply(lpcor, "[", 1))

  #Extract OR and 95% CIs
  library(broom)
  modlist <- lapply(mods, tidy)
  output <- do.call(rbind, modlist) %>% rownames_to_column %>% filter(term == "log2(x)") %>%
    mutate_at(vars(estimate:conf.high), exp)
  output <- if(adj == T){ output %>% mutate(model = "Adj. coffee")
                          } else { 
                          output %>% mutate(model = "Raw") }
}

# Instead generate multiple functions for apply
base <- Caselive_Crs ~ Bmi_C + Alc_Re + Smoke_Stat + alc_drinker_m + L_School + Pa_Mets + 
  Waist_C + strata(Match_Caseset)

mod1 <- function(x) clogit(Caselive_Crs ~ log2(x) + strata(Match_Caseset), data = chcc)
mod2 <- function(x) clogit(update(base, .~. + log2(x)), data = chcc)
mod3 <- function(x) clogit(update(base, .~. + log2(x) + log2(QGE130301+1)), data = chcc)
mod4 <- function(x) clogit(update(base, .~. + log2(x) + Liver_Function_Score), data = chcc)

#Make intensity matrix
ints <- chcc %>% select(Trigonelline:AAMU) %>% as.matrix
fit1 <- apply(ints, 2, mod1)
fit2 <- apply(ints, 2, mod2)
fit3 <- apply(ints, 2, mod3)
fit4 <- apply(ints, 2, mod4)


# Raw model, adjusted covariates, covariates + coffee, covariates + liver fn
CLR.adj <- CLR.coffee(adj = T)
CLR.raw <- CLR.coffee(adj = F)
CLR.cof <- CLR.coffee(adj = T)
CLR.lfn <- CLR.coffee(adj = T)

# Bind output for forest plots: raw vs adjusted, adjusted vs adjusted + coffee
t1a <- bind_rows(CLR.raw, CLR.adj, .id = "Model") %>% arrange(rowname)
t1 <- bind_rows(CLR.adj, CLR.cof) %>% arrange(rowname)

# saveRDS(CLR.all, "OR_hcc_coffee_metabs.rds")
# CLR.all <- readRDS("OR_hcc_coffee_metabs.rds")

# Forest plot of data (raw vs adjusted) ----
rowspace <- rev((1:32)[-seq(3, 32, by = 3)])

#Alternating 18 and 1 point styles
pointvec <- c(rep(c(18, 1), nrow(t1)/2))

dev.off()
# Sets margins (bottom, left, top, right)
par(mar=c(4,4,1,2))

par(mfrow=c(1,2))

library(metafor)
forest(slab = t1$rowname, ilab = t1$model, ilab.xpos = -1.5, ilab.pos = 4,
       x = t1$estimate, ci.ub = t1$conf.high, ci.lb = t1$conf.low, 
       refline = 1, rows = rowspace, 
       xlab = "Odds ratio (per doubling biomarker concentration)",
       ylim = c(1, max(rowspace) + 3), xlim = c(-4, 6), 
       pch = pointvec, psize= rep(1.5, nrow(t1)), cex = 0.8)

text(6, max(rowspace) + 2, "OR [95% CI]", cex = 0.8, pos = 2)
text(c(-4, -1.5), max(rowspace) + 2, c("Compound", "Model"), cex = 0.8, pos=4)

# Other analyses ----------------------------------------------------------------------------------------

# Model coffee only
library(survival)
fit <- clogit(Caselive_Crs ~ QGE130301 + Bmi_C + Alc_Re + Smoke_Stat + alc_drinker_m +
                    L_School + Pa_Mets + Waist_C + strata(Match_Caseset), data = chcc)

# To extract OR and 95% CIs, either use broom or concat numeric values

mod.output <- exp(tidy(caf.CLR)[1, -1])
c(coef(caf.CLR)[1], confint(caf.CLR)[1, ]) %>% exp

library(stargazer)
stargazer(caf.CLR1, type = "text", ci = T, apply.coef = exp)

# CLR using lme4 package
library(lme4)
# gives warning message model nearly unidentifiable
fit.CLR1 <- glmer(Caselive_Crs ~ log2(Caffeine) + Bmi_C + Alc_Re + Smoke_Stat + alc_drinker_m +
            L_School + Pa_Mets + Waist_C + (1 | Match_Caseset), family=binomial, 
            #nAGQ = 17, 
            data = chcc)

summary(fit.CLR)

lapply(mod.list, forest_model, covariates="log2(x)")

# Correlations obesity-coffee biomarkers---------------------------------------------------------

# Join obesity biomarker data and get correlations with coffee biomarkers/intake

library(haven)
ob <- read_sas("sas/obesity_biomarkers.sas7bdat")

chcc <- chcc %>% mutate(Coffee_intake = QGE130301+1)

hcc.ob.case <- chcc %>% select(Idepic, Caselive_Crs, Trigonelline:AAMU, Coffee_intake) %>% 
  left_join(ob, by = "Idepic") %>% filter(Caselive_Crs == 1) %>%
  select(-Idepic, -Idepic_Bio, -Country, -Center)

mat <- as.matrix(hcc.ob[, -1]) %>% log2
mat <- as.matrix(hcc.ob.case[, 1]) %>% log2
mat <- as.matrix(hcc.ob.ctrl[, 1]) %>% log2

cormat <- cor(mat, use = "pairwise.complete.obs")

cof.cormat <- cormat[, "Coffee_intake"]

library(corrplot)
cplot <- corrplot(cormat, method = "square", tl.col = "black", order = "hclust", tl.cex = 0.8,
                  cl.ratio = 0.2, cl.align = "l", cl.pos = "r", mar = c(1,1,1,1))

# Correlations for diketopiperazines
cor.test(mat[, 9], mat[, 19])
cor.test(mat[, 7], mat[, 19])
cor.test(mat[, 8], mat[, 19])
