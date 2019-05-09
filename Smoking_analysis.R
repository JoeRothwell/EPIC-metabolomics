library(dplyr)

meta   <- read.csv("cs_metadata.csv")
allffq <- read.csv("allffq.csv")
ptneg  <- read.delim("EPIC serum RP NEG corrected.txt", skip=4)
ptpos  <- read.csv("EPIC RP POS Corrected final peak table into MPP.csv")

#define metadata and intensities to be used
labs <- meta[meta$present.pos==T, ]
mat  <- t(ptpos[, -c(1:5)])
mat  <- ifelse(mat == 0 | is.na(mat), 1, mat)

df   <- cbind(labs, mat)

filt    <- df %>% filter(smoke == "Yes" | smoke == "No") 
ints.df <- filt %>% select(-(No:fish))
ints    <- as.matrix(ints.df)
ints    <- ifelse(ints == 0 | is.na(ints), 1, ints)
logmat  <- log2(ints)

#filter features present in less than 25% of subjects
filtm   <- logmat[ , colSums(logmat != 0) > 80]

#add smoking 1 or 0 column (may not be necessary)
subj <- filt %>% select(No:fish) %>% mutate(smoker = ifelse(smoke == "Yes", 1, 0))

#logistic regression without metabolomics data
fit <- glm(smoker ~ centre + sex + cof + alcbev + type.plate, data = subj, family=binomial)
summary(fit)

#applying to each spectral feature
smoke.fn <- function(metab) glm(smoker ~ metab + centre + sex + cof + alcbev, data = subj, family=binomial)
log.smoke <- apply(filtm, 2, smoke.fn)

#extract p-values for each feature (ind 2 refers to metab variable, 4 is coefficients line)
pval <- sapply(log.smoke, function(f) summary(f)$coefficients[2, 4])
sort(pval)

