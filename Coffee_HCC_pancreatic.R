library(readr)
library(tidyverse)

# Coffee biomarkers and pancreatic cancer
dat <- read_csv("EPIC pancreas rp pos area or height.csv")
meta <- read_csv("EPIC pancreas case control status.csv")

# Get compound names and subset intensities
cmpds <- dat$`Compound Name`
ints <- t(dat[, -c(1:4)]) #%>% rownames_to_column()
Samples <- rownames(ints)

rownames(ints) <- NULL
df <- cbind(Samples, ints) %>% as_tibble()

# Join metadata to intensities
df1 <- right_join(meta, df, by = "Samples") %>% filter(!is.na(Status)) %>%
  separate(Status_Caseset, into = c("prefix", "match"), sep = "_")
ct <- as.numeric(df1$Status == "Case")

# Subset intensities, scale, plot
mat <- df1[, -(1:5)]
mat <- data.matrix(mat, rownames.force = NA) %>% scale

plot.ts(mat[, 1:8], type = "p")
plot.ts(mat[, 9:15], type = "p")

library(survival)
library(broom)

# Function for CLR and tidy data
multiclr <- function(x) clogit(ct ~ x + strata(match), data = df1)
fits <- apply(mat, 2, multiclr)
results <- map_df(fits, ~tidy(., exponentiate = T)) %>% bind_cols(cmpd = cmpds)

# Categorical analysis
catmat <- apply(mat, 2, function(x) cut_number(x, n = 4, labels = 1:4))
fits1 <- apply(catmat, 2, multiclr)
results1 <- map_df(fits1, ~tidy(., exponentiate = T)) %>% 
  bind_cols(cmpd = rep(cmpds, each = 3))

# Correlation
library(corrplot)
colnames(mat) <- cmpds
cormat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
colnames(cormat) <- NULL
corrplot(cormat, method = "color", order = "hclust")

# HILIC data
dat <- read_csv("EPIC pancreas HILIC pos area.csv")
cmpds <- dat$`Compound Name`
ints <- t(dat[, -c(1:5)]) #%>% rownames_to_column()
Samples <- rownames(ints)

rownames(ints) <- NULL

# Separate identifiers to get joinable sample codes
df <- cbind(Samples, ints) %>% as_tibble() %>% 
  separate(Samples, into = c("prefix", "samp"), sep = "_NR180220_")

meta <- meta %>% separate(Samples, into = c("prefix", "samp"), sep = "_NR171204_")

# Join metadata to intensities
df1 <- inner_join(df, meta, by = "samp", match_fun = str_detect) %>% 
  filter(!is.na(Status)) %>%
  separate(Status_Caseset, into = c("prefix", "match"), sep = "_")
ct <- as.numeric(df1$Status == "Case")

mat <- df1[, 3:7]
mat <- data.matrix(mat, rownames.force = NA) %>% scale

library(survival)
library(broom)

# Function for CLR and tidy data
multiclr <- function(x) clogit(ct ~ x + strata(match), data = df1)
fits <- apply(mat, 2, multiclr)
results <- map_df(fits, ~tidy(., exponentiate = T)) %>% bind_cols(cmpd = cmpds)

# Categorical analysis
catmat <- apply(mat, 2, function(x) cut_number(x, n = 4, labels = 1:4))
fits1 <- apply(catmat, 2, multiclr)
results1 <- map_df(fits1, ~tidy(., exponentiate = T)) %>% 
  bind_cols(cmpd = rep(cmpds, each = 3))

colnames(mat) <- cmpds
cormat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
colnames(cormat) <- NULL
corrplot(cormat, method = "color", order = "hclust")

