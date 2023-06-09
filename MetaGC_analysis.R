# Coffee biomarkers or related compounds and pancreatic cancer in EPIC
# Pancreatic case-control, n=326, 163 cases and controls
library(readr)
library(readxl)
library(tidyverse)
library(haven)
library(amelia)

# Change between Meta-GC folder and VK's writable folder as needed
setwd("/data/Epic_Meta_GC/files")
#setwd("/home/knazev")

# Read in participant data and pos and neg data
dat <- read_sas("meta_gc_caco.sas7bdat")
pos <- read_sas("meta_gc_untg_rp_pos_metabo.sas7bdat")
neg <- read_sas("meta_gc_untg_rp_neg_metabo.sas7bdat")
# Features are in columns

# Check for NAs and lowest values, plot a feature
anyNA(pos) # No NAs
min(sapply(pos$Untg_Rp_Pos_M00001, min)) # Min value in data frame is 15944
plot(pos$Untg_Rp_Pos_M00001)
plot(pos$Untg_Rp_Pos_M00002)

anyNA(neg) # No NAs
min(sapply(neg$Untg_Rp_Neg_M00001, min)) # Min value in data frame is 352917
plot(neg$Untg_Rp_Pos_M00001)
plot(neg$Untg_Rp_Pos_M00002)

# Positive mode----

# Join participant data to metabolomics data
posdat <- right_join(dat, pos, by = "Idepic")

# Subset metabolomics data only using prefix. Log transform data.
posints <- posdat %>% select(starts_with("Untg_Rp_Pos")) %>% as.matrix() %>% log()

# Check for H.pylori status 
table(posdat$HPPOS) # 108 controls and 299 cases

# Function to apply CLR across metabolomics matrix, raw and multivariable-adjusted
# Matching factors are study centre, age at blood collection, sex, fasting status, 
# time of blood collection, menopausal status, exogenous hormone use, phase of menstrual cycle

# time of blood collection, menopausal status, exogenous hormone use, phase of menstrual cycle
# Red meat, processed meat, fruits, vegetables: Qge_0701, Qge_0704, Qge_0401, Qge_02
# Dairy products, fish, eggs, fibre: QgE05 + QgE0801 + QgE0901 + QE_FIBT
# Combined fruits and vegetable variable:
posdat$fv_total <- posdat$QgE0401 + posdat$QgE02

library(survival)
clr.raw <- function(x) clogit(Cncr_Caco_Stom.x ~ x + strata(Match_Caseset.x), data = posdat)

# Model with essential covariates
clr.adj <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + 
                                strata(Match_Caseset.x), data = posdat)

# Essential + red and processed meat, fruit and vegetables
clr.max <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C +
                                QgE0701 + QgE0704 + QgE0401 + QgE02 +
                                strata(Match_Caseset.x), data = posdat)

# Essential + red and processed meat, fruit and vegetables combined, dairy, fish, eggs, fibre
clr.max <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C +
                                QgE0701 + QgE0704 + fv_total +
                                QgE05 + QgE0801 + QgE0901 + QE_FIBT +
                                strata(Match_Caseset.x), data = posdat)

# Model with h.pylori status. First make "unknown" variable so that no participants missing
# this variable are excluded from the model

posdat$HPPOS2 <- as_factor(posdat$HPPOS)
posdat$HPPOS2 <- fct_explicit_na(posdat$HPPOS2)
clr.hpp <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + HPPOS2 +
                                strata(Match_Caseset.x), data = posdat)


# Other potential variables to be adjusted for (from Selenium SAS script ORs - tertiles - BY HPPOS): 
# Red and processed meat, citrus fruit, non-citrus fresh fruit, haem iron intake, fruit and vegetable intake



# Convert categorical variables to factors
# Not sure if samples are matched on fasting status, if so need to remove
varlist <- c("Smoke_Stat", "Alc_Drinker", "Center.x", "L_School", "Fasting_C")
posdat <- posdat %>% mutate(across((varlist), as.factor))


# Run models 
pos.adj <- apply(posints, 2, clr.adj)
pos.max <- apply(posints, 2, clr.max) 
pos.hpp <- apply(posints, 2, clr.hpp)

# Get feature names from data frame label
# Apply attr() across colnames of posints to get vector of feature names
#features <- sapply(pos[ , 1:dim(posints)[2]], attr, "label") %>% unname
features <- colnames(posints)


# Extract raw and adjusted results to data frame and put together
library(broom)
#res.raw <- map_df(pos.raw, tidy, exponentiate = T, conf.int = T) %>% 
#  mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features)

res.adj <- map_df(pos.adj, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(pos.adj))

res.max <- map_df(pos.max, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) #%>% mutate(feat.no = 1:length(pos.adj))

res.hpp <- map_df(pos.hpp, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) #%>% mutate(feat.no = 1:length(pos.adj))

# Make a copy of the results in case of overwrite
res.adj.pos <- res.adj
res.max.pos <- res.max
res.hpp.pos <- res.hpp

# Format OR and CI
tab.adj <- res.adj %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)

tab.max <- res.max %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)

tab.hpp <- res.hpp %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)

# Get FDR adjustment threshold for plot
thr.adj <- res.adj %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max
thr.max <- res.max %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max
thr.hpp <- res.hpp %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max

# Plots
# Fully adjusted model
library(ggrepel)
ggplot(res.adj, aes(x = estimate, y = -log10(p.value))) + geom_point(shape = 1) +
  theme_bw(base_size = 10) +
  xlab("Odds ratio per log increase concentration") + 
  ylab(expression(paste(italic(P), "-value"))) +
  geom_vline(xintercept = 1, linewidth = 0.2, colour = "grey60") + 
  geom_hline(yintercept = -log10(thr.adj), linewidth = 0.2, colour = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 9)) +
  geom_text_repel(aes(label = feature), size = 3, data = res.adj[res.adj$p.value < thr.adj, ]) 

# Fully adjusted model + others
ggplot(res.max, aes(x = estimate, y = -log10(p.value))) + geom_point(shape = 1) +
  theme_bw(base_size = 10) +
  xlab("Odds ratio per log increase concentration") + 
  ylab(expression(paste(italic(P), "-value"))) +
  geom_vline(xintercept = 1, linewidth = 0.2, colour = "grey60") + 
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_hline(yintercept = -log10(thr.max), linewidth = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 9)) +
  geom_text_repel(aes(label = feature), size = 3, data = res.max[res.max$p.value < thr.max, ]) 

# Essential plus H.pylori status adjusted
ggplot(res.hpp, aes(x = estimate, y = -log10(p.value))) + geom_point(shape = 1) +
  theme_bw(base_size = 10) +
  xlab("Odds ratio per log increase concentration") + 
  ylab(expression(paste(italic(P), "-value"))) +
  geom_vline(xintercept = 1, linewidth = 0.2, colour = "grey60") + 
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_hline(yintercept = -log10(thr.hpp), linewidth = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 9)) +
  geom_text_repel(aes(label = feature), size = 3, data = res.hpp[res.hpp$p.value < thr.hpp, ]) 


# Negative mode----

# Join participant data to metabolomics data
negdat <- right_join(dat, neg, by = "Idepic")

# Subset metabolomics data only using prefix. Log transform intensities
negints <- negdat %>% select(starts_with("Untg_Rp_Neg")) %>% as.matrix() %>% log()
dim(negints) # 1130 features

# Function to apply CLR across metabolomics matrix, raw and multivariable-adjusted
# Co-variates to be edited!
library(survival)
#clr.raw <- function(x) clogit(Cncr_Caco_Stom.x ~ x + strata(Match_Caseset.x), data = negdat)

negdat$fv_total <- negdat$QgE0401 + negdat$QgE02

# Model with essential covariates
clr.adj <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + 
                                strata(Match_Caseset.x), data = negdat)

# Essential + red and processed meat, fruit and vegetables
clr.max <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C +
                                QgE0701 + QgE0704 + QgE0401 + QgE02 +
                                strata(Match_Caseset.x), data = negdat)

# Essential + red and processed meat, fruit and vegetables combined, dairy, fish, eggs, fibre
clr.max <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C +
                                QgE0701 + QgE0704 + fv_total +
                                QgE05 + QgE0801 + QgE0901 + QE_FIBT +
                                strata(Match_Caseset.x), data = negdat)

# Model with h.pylori status
negdat$HPPOS2 <- as_factor(negdat$HPPOS)
negdat$HPPOS2 <- fct_explicit_na(negdat$HPPOS2)

clr.hpp <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + HPPOS2 +
                                strata(Match_Caseset.x), data = negdat)

# Small number of missing values in QE_ENERGY (4), QE_ALC (4), L_School (11), Fasting_C (8)

# Convert categorical variables to factors
# Not sure if samples are matched on fasting status, if so need to remove
varlist <- c("Smoke_Stat", "Alc_Drinker", "Center.x", "L_School", "Fasting_C")
negdat <- negdat %>% mutate(across((varlist), as.factor))

# Run models 
neg.adj <- apply(negints, 2, clr.adj)
neg.max <- apply(negints, 2, clr.max) 
neg.hpp <- apply(negints, 2, clr.hpp)

# Get feature names from data frame label
# Apply attr() across colnames of posints to get vector of feature names
features <- sapply(neg[ , 1:dim(negints)[2]], attr, "label") %>% unname
features <- colnames(negints)

# Extract raw and adjusted results to data frame and put together
library(broom)
#res.raw <- map_df(neg.raw, tidy, exponentiate = T, conf.int = T) %>% 
#  mutate(p.adj = p.adjust(p.value, method = "fdr"))

res.adj <- map_df(neg.adj, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(neg.adj))

res.max <- map_df(neg.max, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term)

res.hpp <- map_df(neg.hpp, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term)

# Make a copy of the results in case of overwrite
res.adj.neg <- res.adj


# Format OR and CI
tab.adj <- res.adj %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)

tab.max <- res.max %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)

tab.hpp <- res.hpp %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)

# Get FDR adjustment threshold for plot
thr.adj <- res.adj %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max
thr.max <- res.max %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max
thr.hpp <- res.hpp %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max

# Plots
# Fully adjusted model
library(ggrepel)
ggplot(res.adj, aes(x = estimate, y = -log10(p.value))) + geom_point(shape = 1) +
  theme_bw(base_size = 10) +
  xlab("Odds ratio per log increase concentration") + 
  ylab(expression(paste(italic(P), "-value"))) +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = -log10(thr.adj), size = 0.2, colour = "grey60") +
  geom_hline(yintercept = 0, size = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 9)) +
  geom_text_repel(aes(label = feature), size = 3, data = res.adj[res.adj$p.value < thr.adj, ]) 

# Fully adjusted model + others
ggplot(res.max, aes(x = estimate, y = -log10(p.value))) + geom_point(shape = 1) +
  theme_bw(base_size = 10) +
  xlab("Odds ratio per log increase concentration") + 
  ylab(expression(paste(italic(P), "-value"))) +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = 0, size = 0.2, colour = "grey60") +
  geom_hline(yintercept = -log10(thr.max), size = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 9)) +
  geom_text_repel(aes(label = feature), size = 3, data = res.max[res.max$p.value < thr.max, ]) 

# Fully adjusted model + others
ggplot(res.hpp, aes(x = estimate, y = -log10(p.value))) + geom_point(shape = 1) +
  theme_bw(base_size = 10) +
  xlab("Odds ratio per log increase concentration") + 
  ylab(expression(paste(italic(P), "-value"))) +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = 0, size = 0.2, colour = "grey60") +
  geom_hline(yintercept = -log10(thr.hpp), size = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 9)) +
  geom_text_repel(aes(label = feature), size = 3, data = res.hpp[res.hpp$p.value < thr.hpp, ])


# Correlations----
# First rename res.adj to res.adj.pos and res.adj.neg as appropriate
sig.posfeats <- res.adj.pos %>% filter(p.adj < 0.05) %>% pull(feat.no)
sig.negfeats <- res.adj.neg %>% filter(p.adj < 0.05) %>% pull(feat.no)

# Get intensities of discriminants from peak tables and bind together
posints.sig <- posints[, sig.posfeats]
negints.sig <- negints[, sig.negfeats]
logmat <- cbind(posints.sig, negints.sig) %>% log10

# Get feature names
posnames <- paste("Neg", res.adj.pos[sig.posfeats, ]$feature)
negnames <- paste("Pos", res.adj.neg[sig.negfeats, ]$feature)
sigfeatnames <- c(posnames, negnames)


library(corrplot)
cormat <- cor(logmat, method = "pearson")
colnames(cormat) <- sigfeatnames
rownames(cormat) <- sigfeatnames
corrplot(cormat, method = "square", order = "hclust", tl.col = "black")

# Get cluster order (assign corrplot to corp first)
clust.ord <- unique(corp$corrPos$xName)
cormat1 <- cormat[clust.ord, clust.ord]

write.csv(cormat1, "metabolite correlations MetaGC.csv")
