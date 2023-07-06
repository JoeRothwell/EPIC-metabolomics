# MetaGC sensitivity analyses (as suggested by Mazda). 
# Run appropriate commands first from MetaGC_hp_subsets.R to get each relevant data subset

# (a)	run sensitivity analyses only on the n=105 case-control pairs with Hppos = 1,
# Subset metabolomics data only using prefix. Log transform data.
# Model with essential covariates
posints <- concordant11 %>% select(starts_with("Untg_Rp_Pos")) %>% as.matrix() %>% log()

clr.adj <- function(x) clogit(Cncr_Caco_Stom.x ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + 
                                strata(Match_Caseset.x), data = concordant11)

library(survival)
varlist <- c("Smoke_Stat", "Alc_Drinker", "Center.x", "L_School", "Fasting_C")
concordant11 <- concordant11 %>% mutate(across((varlist), as.factor))
# Apply model across subset
pos.adj <- apply(posints, 2, clr.adj)

features <- colnames(posints)

library(broom)
res.adj <- map_df(pos.adj, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(pos.adj))



# (b)	run sensitivity analyses on the n=105 + n=71 (i.e., n=176) case-control pairs where the case is Hppos =1, 
# but the control could be 1 or 0.

# This is concordant11 + discordant10
discordant10 <- discordant10 %>% select( -disc.type)
caseposall <- bind_rows(concordant11, discordant10)


# (c)	Possibly consider an analysis on the n=11 + n=17 where the cases are Hppos negative, but this will be weak.


# (d)	Break matching and run unconditional analyses by Hppos, but with adjustment for country instead of by centre.


# (e)	In all sub-group analyses a, b, c and d above, also run analyses on the subset of n=230 who have missing data, 
# possibly with heterogeneity analyses?
  