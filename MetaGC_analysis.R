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

# Read in participant data and pos and neg data, removing unneeded cols from feature tables
dat <- read_sas("meta_gc_caco.sas7bdat")
pos <- read_sas("meta_gc_untg_rp_pos_metabo.sas7bdat") %>% select(-(Idepic_Bio:Cncr_Caco_Stom))
neg <- read_sas("meta_gc_untg_rp_neg_metabo.sas7bdat") %>% 
  select(-(Idepic_Bio:Center)) %>% select(-Match_Caseset, -Cncr_Caco_Stom, -X)
# Features are in columns

# Concatenate pos and neg features
posneg <- bind_cols(pos, neg)

# Check for NAs and lowest values, plot a feature
anyNA(posneg) # No NAs
min(sapply(posneg$Untg_Rp_Neg_M00001, min)) # Min value in data frame is 352917
plot(posneg$Untg_Rp_Pos_M00001)
plot(posneg$Untg_Rp_Pos_M00002)

# Join participant data to metabolomics data (note, data still called posdat though neg features included)
posdat <- right_join(dat, posneg, by = "Idepic")

# Subset metabolomics data only using prefix. Log2 transform data.
posints <- posdat %>% select(starts_with("Untg_Rp")) %>% as.matrix() %>% log2()

# Get feature names from data frame label
# Apply attr() across colnames of posints to get vector of feature names
#features <- sapply(pos[ , 1:dim(posints)[2]], attr, "label") %>% unname
features <- colnames(posints)

# Function to apply CLR across metabolomics matrix. Matching factors are study centre, sex, age at blood 
# collection, time of blood collection, fasting status, menopausal status, exogenous hormone use, phase 
# of menstrual cycle

# Convert categorical variables to factors
varlist <- c("Smoke_Stat", "Alc_Drinker", "Center", "L_School", "Fasting_C")
posdat <- posdat %>% mutate(across((varlist), as.factor))

# Foods to potentially adjust for:
# Red meat, processed meat, fruits, vegetables: Qge_0701, Qge_0704, Qge_0401, Qge_02
# Dairy products, fish, eggs, fibre: QgE05 + QgE0801 + QgE0901 + QE_FIBT
# Combined fruits and vegetable variable:
posdat$fv_total <- posdat$QgE0401 + posdat$QgE02

# Other potential variables to be adjusted for (from Selenium SAS script ORs - tertiles - BY HPPOS): 
# Red and processed meat, citrus fruit, non-citrus fresh fruit, haem iron intake, fruit and vegetable intake

library(survival)
# Raw model (not used)
#clr.raw <- function(x) clogit(Cncr_Caco_Stom ~ x + strata(Match_Caseset.x), data = posdat)

#res.raw <- map_df(pos.raw, tidy, exponentiate = T, conf.int = T) %>% 
#  mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features)

# Define model with essential covariates
clr.adj <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + 
                                strata(Match_Caseset), data = posdat)

# Run models, extract raw and adjusted results to data frame and put together
mod.adj <- apply(posints, 2, clr.adj)

library(broom)
res.adj <- map_df(mod.adj, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.adj))

# Format OR and CI
tab.adj <- res.adj %>%
  select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)



# Model 2: Essential + red and processed meat, fruit and vegetables (to be replaced by mod3)
clr.max <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C +
                                QgE0701 + QgE0704 + QgE0401 + QgE02 +
                                strata(Match_Caseset), data = posdat)

# Apply across matrix 
mod.max <- apply(posints, 2, clr.max) 

# Extract results
res.max <- map_df(mod.max, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.max))

# Format OR and CI
tab.max <- res.max %>%
  select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)



# Model 3: Essential + red and processed meat, fruit and vegetables combined, dairy, fish, eggs, fibre
clr.max <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C +
                                QgE0701 + QgE0704 + fv_total +
                                QgE05 + QgE0801 + QgE0901 + QE_FIBT +
                                strata(Match_Caseset), data = posdat)

# Apply across matrix 
mod.max <- apply(posints, 2, clr.max) 

# Extract results
res.max <- map_df(mod.max, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.max))

# Format OR and CI
tab.max <- res.max %>%
  select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)


# Model with h.pylori status. First make "unknown" variable so that no participants missing
# this variable are excluded from the model

# Check for H.pylori status 
table(posdat$HPPOS) # 108 controls and 299 cases

posdat$HPPOS2 <- as_factor(posdat$HPPOS)
posdat$HPPOS2 <- fct_explicit_na(posdat$HPPOS2)
clr.hpp <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + HPPOS2 +
                                strata(Match_Caseset), data = posdat)

# Run models 
mod.hpp <- apply(posints, 2, clr.hpp)

res.hpp <- map_df(mod.hpp, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.hpp))

tab.hpp <- res.hpp %>%
  select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, conf.low, hyph, conf.high, B2, sep = "") %>% filter(p.adj < 0.05)


# Plots. First need to calculate FDR p-value thresholds
thr.adj <- res.adj %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max
thr.max <- res.max %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max
thr.hpp <- res.hpp %>% filter(p.adj < 0.05) %>% pull(p.value) %>% max


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



# Correlations----
# *To be updated for pos and neg togther*
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
