# Coffee biomarkers or related compounds and pancreatic cancer in EPIC
# Pancreatic case-control, n=326, 163 cases and controls
library(readr)
library(readxl)
library(tidyverse)
library(haven)

# Coffee biomarkers and pancreatic cancer
dat <- read_csv("EPIC pancreas rp pos area or height.csv")

# Case-control status
meta <- read_csv("EPIC pancreas case control status.csv") 
# ID conversion data from Pekka with Idepic_Bio
ids <- read_xlsx("EPIC pancreas case contraol status 160221.xlsx")

# Get subject data (from Bertrand Feb 2021)
subjects <- read_tsv("Panc_CaCo/panc_caco_uc.txt") %>%
  mutate_at(vars(L_School, Smoke_Stat, Pa_Total, Fasting_C), as.factor)

ids.subj <- left_join(ids, subjects, by = "Idepic_Bio")

# Get compound names and subset intensities
cmpds <- paste(dat$`Compound Name`, "_int", sep = "")
ints <- t(dat[, -c(1:4)]) 
colnames(ints) <- cmpds
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

plot.ts(mat[, 1:8], type = "p", main = "AAMU - Glycocholic acid")
plot.ts(mat[, 9:15], type = "p", main = "Hypoxanthine - Cyclo(iso-pro)")

library(survival)
library(broom)

# Function for CLR and tidy data. Base model
multiclr <- function(x) clogit(ct ~ x + strata(match), data = df1)

# Adjusted model
multiclr <- function(x) clogit(ct ~ x + Bmi_C + Pa_Total + Smoke_Stat + QE_ENERGY + Fasting_C +
                                  QE_ALC + L_School + strata(match), data = df1)

fits <- apply(mat, 2, multiclr)
results <- map_df(fits, ~tidy(., exponentiate = T, conf.int = T)) %>% filter(term == "x") %>%
  bind_cols(cmpd = cmpds)

RR.ci(data = results)


# Categorical analysis
catmat <- apply(mat, 2, function(x) cut_number(x, n = 4, labels = 1:4))
fits1 <- apply(catmat, 2, multiclr)

# Tidy results into data frame. Unadjusted model:
results1 <- map_df(fits1, ~tidy(., exponentiate = T, conf.int = T)) %>% 
  bind_cols(cmpd = rep(cmpds, each = 3))

# Adjusted model:
results1 <- map_df(fits1, ~tidy(., exponentiate = T, conf.int = T)) %>% 
  filter(str_detect(term, "x")) %>%
  bind_cols(cmpd = rep(cmpds, each = 3))

# Format
tabcat <- results1 %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "")

# Correlation
library(corrplot)
colnames(mat) <- cmpds
cormat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
colnames(cormat) <- NULL
corrplot(cormat, method = "color", order = "hclust")


### HILIC data
dat <- read_csv("EPIC pancreas HILIC pos area.csv")

# Case-control status
meta <- read_csv("EPIC pancreas case control status.csv") 
meta <- meta %>% separate(Samples, into = c("prefix", "samp"), sep = "_NR171204_")

cmpds <- paste(dat$`Compound Name`, "_H_int", sep = "")
ints <- t(dat[, -c(1:5)]) #%>% rownames_to_column()
Samples <- rownames(ints)

rownames(ints) <- NULL

# Separate identifiers to get joinable sample codes
df <- cbind(Samples, ints) %>% as_tibble() %>% 
  separate(Samples, into = c("prefix", "samp"), sep = "_NR180220_", remove = F)

#meta <- meta %>% separate(Samples, into = c("prefix", "samp"), sep = "_NR171204_")

# Join metadata to intensities
df1 <- inner_join(df, meta, by = "samp", match_fun = str_detect) %>% 
  filter(!is.na(Status)) %>%
  separate(Status_Caseset, into = c("prefix", "match"), sep = "_") %>%
  # Use str_replace to get IDs that are joinable to the participant data
  mutate(Samples1 = str_replace(Samples, "NR180220", "NR171204"))

# Vector of case-control status
ct <- as.numeric(df1$Status == "Case")

# Join participant data
df2 <- inner_join(df1, ids.subj, by = c("Samples1"="Comments"))

mat <- df1[, 4:8]
mat <- data.matrix(mat, rownames.force = NA) %>% scale

library(survival)
library(broom)

# Function for CLR and tidy data
multiclr <- function(x) clogit(ct ~ x + strata(match), data = df1)

multiclr <- function(x) clogit(ct ~ x + Bmi_C + Pa_Total + Smoke_Stat + QE_ENERGY + Fasting_C +
                                 QE_ALC + L_School + strata(match), data = df2)

fits <- apply(mat, 2, multiclr)
results <- map_df(fits, ~tidy(., exponentiate = T, conf.int = T)) %>% 
  filter(str_detect(term, "x")) %>% bind_cols(cmpd = cmpds)

# Format
tabcat <- results %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "")

# Categorical analysis
catmat <- apply(mat, 2, function(x) cut_number(x, n = 4, labels = 1:4))
fits1 <- apply(catmat, 2, multiclr)

# Extract results to data frame
results1 <- map_df(fits1, ~tidy(., exponentiate = T, conf.int = T)) %>% 
  filter(str_detect(term, "x")) %>% bind_cols(cmpd = rep(cmpds, each = 3))

# Format
tabcat <- results1 %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("OR", estimate, B1, conf.low, hyph, conf.high, B2, sep = "")


colnames(mat) <- cmpds
cormat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
colnames(cormat) <- NULL
corrplot(cormat, method = "color", order = "hclust")

