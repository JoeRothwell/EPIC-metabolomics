# Coffee biomarkers or related compounds and pancreatic cancer in EPIC
# Pancreatic case-control, n=326, 163 cases and controls
library(readr)
library(readxl)
library(tidyverse)
library(haven)

# Read in intensity data and harmonise colnames between RP and HILIC
dat <- read_csv("EPIC pancreas rp pos area or height.csv")
dat1 <- read_csv("EPIC pancreas HILIC pos area.csv") %>%
  mutate(`Compound Name` = paste(`Compound Name`, "_H", sep = ""))

# Replace part of ID to give samples same name in RP and HILIC
colnames(dat1) <- str_replace(colnames(dat1), "NR180220", "NR171204")
rp    <- dat[, intersect(colnames(dat), colnames(dat1))]
hilic <- dat1[, intersect(colnames(dat), colnames(dat1))]
# Bind together RP and HILIC data
dat <- bind_rows(rp, hilic)

# Case-control status.
meta <- read_csv("EPIC pancreas case control status.csv")
# ID conversion data from Pekka with Idepic_Bio
ids <- read_xlsx("EPIC pancreas case contraol status 160221.xlsx")

# Get subject data (from Bertrand Feb 2021)
subjects <- read_tsv("Panc_CaCo/panc_caco_uc.txt") %>%
  mutate_at(vars(L_School, Smoke_Stat, Pa_Total, Fasting_C), as.factor)

ids.subj <- left_join(ids, subjects, by = "Idepic_Bio")

# Subset intensities and add compound names
ints <- dat %>% select(starts_with("EPIC_")) %>% t
#ints <- t(dat[, -c(1:4)]) 
cmpds <- paste(dat$`Compound Name`, "_int", sep = "")
colnames(ints) <- cmpds

# Convert the sample names to a df column
Samples <- rownames(ints)
rownames(ints) <- NULL
df <- cbind(Samples, ints) %>% as_tibble()

# Join case-control status to intensities
df1 <- inner_join(meta, df, by = "Samples") %>% filter(!is.na(Status)) %>%
  separate(Status_Caseset, into = c("prefix", "match"), sep = "_")

# Join to IDs and participant data
df2 <- inner_join(df1, ids.subj, by = c("Samples"="Comments"))
df1 <- df2
  
# Get vector of case-control status  
ct <- as.numeric(df1$Status == "Case")

# Subset intensities, scale, plot
mat <- df1 %>% select(ends_with("_int"))
mat <- data.matrix(mat, rownames.force = NA) %>% scale
colnames(mat) <- abbreviate(cmpds, use.classes = F)

#plot.ts(mat[, 1:8], type = "p", main = "AAMU - Glycocholic acid")
#plot.ts(mat[, 9:15], type = "p", main = "Hypoxanthine - Cyclo(iso-pro)")


# Run conditional logistic models for 20 biomarkers
library(survival)
library(broom)

# Function for CLR and tidy data. Base model
multiclr1 <- function(x) clogit(ct ~ x + strata(match), data = df1)

# Adjusted model
multiclr2 <- function(x) clogit(ct ~ x + Bmi_C + Pa_Total + Smoke_Stat + QE_ENERGY + 
                                 Fasting_C + QE_ALC + L_School + strata(match), data = df1)
# Run models
fits2 <- apply(mat, 2, multiclr2)
fits1 <- apply(mat, 2, multiclr1)

# Extract raw and adjusted results to data frame and put together
results <- map_df(c(fits1, fits2), ~tidy(., exponentiate = T, conf.int = T)) %>% 
  filter(term == "x") %>% bind_cols(mod = c(rep("raw", 20), rep("adjusted", 20)), 
                                            cmpd = rep(cmpds, 2))

# Format OR and CI
tabcat <- results %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "")

# Categorical analysis
catmat <- apply(mat, 2, function(x) cut_number(x, n = 4, labels = 1:4))
fits3 <- apply(catmat, 2, multiclr1)
fits4 <- apply(catmat, 2, multiclr2)

# Extract raw and adjusted results to data frame and put together
results1 <- map_df(c(fits3, fits4), ~tidy(., exponentiate = T, conf.int = T)) %>% 
  filter(str_detect(term, "x")) %>% 
  mutate(mod = rep(c("raw", "adjusted"), each = 60), cmpd = rep(rep(cmpds, each = 3), 2))

# Format OR and CI
tabcat <- results1 %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "")

# Correlations
colnames(mat) <- cmpds
cormat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
colnames(cormat) <- NULL
corrplot(cormat, method = "color", order = "hclust")

