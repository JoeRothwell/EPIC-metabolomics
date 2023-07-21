# MetaGC sensitivity analyses (as suggested by Mazda). 
# Run appropriate commands first from MetaGC_hp_subsets.R to get each relevant data subset

# (a)	run sensitivity analyses only on the n=105 case-control pairs with Hppos = 1,
# Subset metabolomics data only using prefix. Log transform data.
# Model with essential covariates
library(tidyverse)
posints <- concordant11 %>% select(starts_with("Untg_Rp_Pos")) %>% as.matrix() %>% log2()

clr.adj <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + QE_ENERGY + QE_ALC + 
                                L_School + Fasting_C + strata(Match_Caseset), data = concordant11)

library(survival)
varlist <- c("Smoke_Stat", "Alc_Drinker", "Center", "L_School", "Fasting_C")
concordant11 <- concordant11 %>% mutate(across((varlist), as.factor))
# Apply model across subset
mod.adj2 <- apply(posints, 2, clr.adj)

features <- colnames(posints)

library(broom)
res.adj2 <- map_df(mod.adj2, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.adj2))

# No significant features after p-value adjustment




# (b)	run sensitivity analyses on the n=105 + n=71 (i.e., n=176) case-control pairs where the case is Hppos=1, 
# but the control could be 1 or 0.

# This is concordant11 + discordant10 (168 + 134 obs = actually 151 pairs)
discordant10 <- discordant10 %>% mutate(across((varlist), as.factor)) %>% select(-disc.type)
caseposall <- bind_rows(concordant11, discordant10)
posints <- caseposall %>% select(starts_with("Untg_Rp_Pos")) %>% as.matrix() %>% log2()

mod.adj3 <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + 
                                strata(Match_Caseset), data = caseposall)

# Apply model across subset
mod.adj3 <- apply(posints, 2, res.adj3)

res.adj3 <- map_df(mod.adj3, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.adj3))
# No significant features after p-value adjustment




# (c)	Possibly consider an analysis on the n=11 + n=17 where the cases are Hppos negative, but this will be weak.

# This is concordant00 + discordant01 (28 pairs)
casenegall <- bind_rows(concordant00, discordant01)
posints <- casenegall %>% select(starts_with("Untg_Rp_Pos")) %>% as.matrix() %>% log2()

mod.adj4 <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + 
                                strata(Match_Caseset), data = casenegall)

# Apply model across subset
pos.adj <- apply(posints, 2, mod.adj4)

res.adj <- map_df(pos.adj, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(pos.adj))




# (d)	Break matching and run unconditional analyses by Hppos, but with adjustment for country instead of by centre.

# Matching factors are study centre, age at blood collection, sex, fasting status, 
# time of blood collection, menopausal status, exogenous hormone use, phase of menstrual cycle
# Model with essential covariates

# Better to leave unknown HPPOS status
#posdat$HPPOS2 <- as_factor(posdat$HPPOS)
#posdat$HPPOS2 <- fct_explicit_na(posdat$HPPOS2)

# Unmatched LR model adjusting additionally for Hppos and matching factors where possible: sex, country,
# hormone use, fasting status
glm.hpp <- function(x) glm(Cncr_Caco_Stom ~ x + Sex + Country + #Use_Horm + 
                             Bmi_C + Smoke_Stat + Pa_Total + 
                                Age_Blood + QE_ENERGY + QE_ALC + L_School + Fasting_C + HPPOS,
                                data = posdat, family = "binomial")

# Reset posdat from MetaGC_hp_subsets.R
posints <- posdat %>% select(starts_with("Untg_Rp")) %>% as.matrix() %>% log2()

features <- colnames(posints)

# Convert categorical variables to factors
varlist <- c("Smoke_Stat", "Alc_Drinker", "Country", "L_School", "Fasting_C", "Use_Horm", "Sex", "HPPOS")
posdat <- posdat %>% mutate(across((varlist), as.factor))

# Apply model across subset
pos.adj <- apply(posints, 2, glm.hpp)

#res.adj <- map_df(pos.adj, tidy, exponentiate = T, conf.int = T) %>% 
res.adj <- map_df(pos.adj, tidy) %>% 
  filter(term == "x") %>% 
  mutate(OR = exp(estimate), p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(pos.adj))



# (e)	In all sub-group analyses a, b, c and d above, also run analyses on the subset of n=230 who have missing data, 
# possibly with heterogeneity analyses? [to do]
  