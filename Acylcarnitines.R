library(tidyverse)

# Filter compounds by a threshold of detection
ac <- read_csv("acylcarnitines.csv") %>% 
  filter(`Detection Frequency` > 50)

# Get acylcarnitine names from metadata
acnames <- ac$Name
acCno   <- ac$`C number`

# Prepare matrix of intensities and generate correlation matrix. Remove metadata and sample repeats
ints <- ac %>% select(contains("[Area]"), -contains("repeat"))
colnames(ints)

mat              <- as.matrix(ints) %>% t
colnames(mat)    <- acnames
mat[mat == 0]    <- NA
cormat           <- cor(mat, use = "pairwise.complete.obs")
rownames(cormat) <- acnames
colnames(cormat) <- NULL

# Correlations between individual ACs
library(corrplot)
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.6)

# Correlations with foods. Join FFQ data and WCRF scores to subject metadata
ffq <- read_csv("allffq.csv")
meta  <- read_csv("cs_metadata.csv")
library(haven)
wcrf <- read_dta("data/Wcrf_Score.dta") %>% select(Idepic, Wcrf_C_Cal)
meta.ffq <- meta %>% left_join(ffq, by="Idepic")
meta.ffq <- meta.ffq %>% left_join(wcrf, by = "Idepic")

#Extract the data file numbers from the column names for joining
sample.list <- colnames(ints)
datalabel <- sample.list %>% str_match_all("[0-9]+") %>% unlist %>% as.numeric


#Join AC intensities and FFQ/metadata
int.df <- cbind(datalabel = paste(datalabel, "(raw)", sep=""), mat) %>% tbl_df
all.df <- left_join(int.df, meta.ffq, by = "datalabel") %>%
   filter(Wcrf_C_Cal != 3) %>% mutate(score_cat = ifelse(Wcrf_C_Cal %in% 1:2, 0, 1))

#Convert data frame AC subset to numeric values for correlation and model
acnmat <- select(all.df, contains(":")) %>% data.matrix
logmat <- log2(acnmat)
ffqmat <- select(all.df, starts_with("QGE"))
meta1 <- select(all.df, order:score_cat)

#test for missing values
library(Amelia)
missmap(all.df[, 142:381])
colSums(ffqmat > 0)

#Exclude zero SD columns
nonzerosd <- which(apply(ffqmat, 2, function(x) sd(x, na.rm=T) != 0))
#Exclude foods recorded in less than 10 subjects
filt <- colSums(ffqmat > 0) >10
ffqmat <- ffqmat[, filt]

accor <- cor(acnmat, ffqmat, use = "pairwise.complete.obs")
#NAs in QGE06030202 (176)

library(gplots)
heatmap.2(accor, trace = "none", col=redblue(256), Colv = F)
max(accor)
#get highest correlation
which(accor == max(accor), arr.ind = TRUE)


# Association with WCRF scores ----

varlist <- c("type.plate", "country", "sex")
meta1 <- meta1 %>% mutate_at(vars(varlist), as.factor)

meta1
#logistic regression scores 1 and 2 vs 4 and 5
glm.ac <- function(x) glm(score_cat ~ x + type.plate + country + sex, data = meta1, family = "binomial")

# apply function across metabolite matrix
multifit <- apply(logmat, 2, glm.ac)

library(broom)
p <- map_df(multifit, tidy) %>% filter(term == "x") %>% mutate(Cmpd = acnames, Cnumber = acCno)
p <- mutate(p, p.fdr = p.adjust(p.value, method = "fdr"))

# calculation of fold changes for volcano plot. Make sure to use non-log data
df       <- data.frame(scorecat = meta1$score_cat, acnmat)
means    <- aggregate(. ~ scorecat, data = df, mean) %>% t
# Get fold change of high score over low score
meanfc   <- means[, 2]/means[, 1]

df <- data.frame(p, Fold_change = meanfc[-1]) %>% 
  mutate(direction   = ifelse(Fold_change > 1, "high_score", "low_score"))

library(ggplot2)
ggplot(df, aes(x = reorder(Cmpd, -p.fdr), y = -log10(p.value), shape=direction, colour=direction)) + 
  theme_minimal(base_size = 10) +
  geom_point(show.legend = F) + 
  geom_hline(yintercept = 3, linetype = "dashed") +
  ylab("-log10(FDR-adjusted p-value)") + xlab("Compound") +
  facet_grid(. ~ Cnumber, scales = "free_x", space = "free_x", switch= "y") +
  theme(axis.text.x = element_text(angle = 90, size=7, hjust = 0.95, vjust = 0.5))
#ggtitle("Metabolite associations with WCRF score (cal)") +
#ggsave("WCRF score associations FAs.svg")





