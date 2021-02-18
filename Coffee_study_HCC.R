library(readr)
library(tidyverse)
library(haven)

# Coffee biomarkers and pancreatic and liver cancer 
meta <- read_sas("merged_untarg1.sas7bdat") %>% 
  select(Id_Bma, Caselive_Crs, Match_Caseset, Bmi_C, L_School, Pa_Mets_C, QE_ALC, Smoke_Stat,
         Fasting_C, Waist_C, Alc_Re, Alc_Drinker, Idepic) %>%
  mutate_at(vars(L_School, Pa_Mets_C, Smoke_Stat, Fasting_C, Alc_Drinker), as.factor) %>%
  mutate(Id_Bma = str_sub(Id_Bma, start = 10L, end = 19L))

lf <- read_sas("live_caco1.sas7bdat") %>% 
  select(Idepic, Liver_Function_Score, Liver_Function_Score_C)

meta <- left_join(meta, lf, by = "Idepic")

# Compound data (add _H for hilic)
dat1 <- read_csv("EPIC HCC rp pos area or height.csv")
dat2 <- read_csv("EPIC HCC HILIC pos area.csv") %>% 
  mutate(`Compound Name` = str_c(`Compound Name`, "_H"))

# Standardise colnames to bind rows together
colnames(dat2) <- NULL
colnames(dat2) <- colnames(dat1) 
dat <- rbind(dat1, dat2)

# Get compound names and subset intensities
cmpds <- paste(dat$`Compound Name`, "_int", sep = "")
ints <- t(dat[, -c(1:4)]) #%>% rownames_to_column()
Id_Bma <- rownames(ints)

rownames(ints) <- NULL
colnames(ints) <- cmpds
df <- cbind(Id_Bma, ints) %>% as_tibble() %>% 
  mutate(Id_Bma = str_sub(Id_Bma, start = -14L, end = -5L))

# Join metadata to intensities
df1 <- inner_join(meta, df, by = "Id_Bma")
#ct <- as.numeric(df1$Status == "Case")

# Subset intensities, scale, plot
mat <- df1 %>% select(ends_with("_int"))
#mat <- df1[, -(1:3)]
mat <- data.matrix(mat, rownames.force = NA) %>% scale

plot.ts(mat[, 1:10], type = "p")
plot.ts(mat[, 11:20], type = "p")

library(survival)
library(broom)

# Function for CLR and tidy data
multiclr <- function(x) clogit(Caselive_Crs ~ x + strata(Match_Caseset), data = df1)

# Adjusted model
multiclr <- function(x) clogit(Caselive_Crs ~ x + Bmi_C + L_School + #Pa_Mets_C + 
            QE_ALC + Smoke_Stat + Fasting_C + strata(Match_Caseset), data = df1)

# Adjusted model + liver function score (doesn't work)
multiclr <- function(x) clogit(Caselive_Crs ~ x + Bmi_C + L_School + #Pa_Mets_C + 
            QE_ALC + Smoke_Stat + Fasting_C + Liver_Function_Score + strata(Match_Caseset), 
            data = df1)

fits <- apply(mat, 2, multiclr)
results <- map_df(fits, ~tidy(., exponentiate = T, conf.int = T)) %>% filter(term == "x") %>%
  bind_cols(cmpd = cmpds)
  
RR.ci(data = results)

# Categorical analysis
catmat <- apply(mat, 2, function(x) cut_number(x, n = 4, labels = 1:4))
fits1 <- apply(catmat, 2, multiclr)
results1 <- map_df(fits1, ~tidy(., exponentiate = T)) %>% 
  bind_cols(cmpd = rep(cmpds, each = 3))
