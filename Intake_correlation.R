# Finds features associated with food intake by partial correlation
# Coffee QGE130301; Red meat QGE0701; Fish QGE0801; Fruits, nuts and seeds QGE04; Alc beverages QGE14;
# Beer 140301; Total dietary fibre QEFIBT

### Function for baseline correction -------------------------------------------------------------------
baselines <- function(x) {
  
  #calculation of raw and corrected peaks for new data 4-4-2017
  #after raw data check, pos mode@ 357 is classified as raw, 204 is corr
  
  #use old data colnames to determine baseline heights (high/low) for each sample
  # get colnames from feature tables
  
  library(tidyverse)
  posname <- read_csv("Old/EPIC RP POS Corrected final peak table into MPP.csv") %>% colnames
  negname <- read_tsv("Old/EPIC serum RP NEG corrected.txt", skip=4) %>% colnames
  
  library(stringr)
  # get neg sample numbers and their baseline levels #data_frame doesn't convert to factor
  sampno <- negname %>% str_match_all("[0-9]+") %>% unlist %>% as.numeric
  negbl <- negname %>% str_match("corr|raw")
  negdf  <- data_frame(sampno, negbl = negbl[-1])
  
  #data frame of positive baseline levels
  sampno <- posname %>% str_match_all("[0-9]+") %>% unlist %>% as.numeric
  posbl <- posname %>% str_match("corr|raw")
  posdf  <- data_frame(sampno, posbl = posbl[-(1:5)])
  
  bl <- merge(negdf, posdf, by="sampno", all.x = T)
  
  which(is.na(bl$posbl))
  bl$posbl[c(204, 355)] <- c("corr", "raw")
  
  #insert NAs into baseline cols and subset
  bl <- read_csv("cs_metadata.csv") %>% bind_cols(bl) %>% 
    mutate(pos.bl = if_else(stype == "SA", posbl, country)) %>%
    mutate(neg.bl = if_else(stype == "SA", negbl, country)) %>%
    select(pos.bl, neg.bl)
  #return(bl)
  # saveRDS(bl, "baselines pos neg.rds")
}
bl <- baselines()

### Function for partial correlation analysis of metabolites ------------------------------------------

intake.corr <- function(food, pos = T, incr = T, impute = F, min.sample = 340, pcutoff = 0.05) {

  require(tidyverse)
  #see baseline correction.R for details of baseline co-variate
  baselines <- readRDS("baselines pos neg.rds")

  ### Define and process sample metadata (food, lifestyle, technical) ###
  #Get alcohol (g) and BMI data
  alc_g <- read.csv("alcohol/alcohol.csv") %>% select(Idepic, Qe_Alc, R_BMI)

  # Read in main metadata, add alcohol and BMI data and log transformed intakes
  meta  <- read.csv("cs_metadata.csv") %>% left_join(alc_g, by="Idepic") %>%
    mutate(logcof = log(cof + 1), logalc = log(alcbev + 1), logredmeat = log(redmeat + 1),
         logprocmeat = log(procmeat + 1), logfish = log(fish + 1), logQe_Alc = log(Qe_Alc + 1))

  #add cup volumes in ml: France, 146.59; Italy, 55.2; Greece, 135.48; Germany, 209.32
  cupvols <- data.frame(country = levels(meta$country), cupvol = c(146.59, 209.32, 135.48, 55.2))
  meta    <- meta %>% left_join(cupvols, by="country") %>% mutate(cups = cof/cupvol)

  # Filter the 498 obs to get the 451 measured subjects only (excluding those not properly injected in pos)
  # Excluded samples are 204 and 355
  sampvec <- meta$present.pos == T & meta$stype == "SA"

  ### Read and process metabolomics data ------------------------------------------------------------
  
  if(pos == T) {  
    ion.mode <- "Pos"
    meta$baseline <- baselines$pos.bl
    pt <- read.delim("data/EPIC Cross sectional RP POS Feature table.txt", skip=4, row.names = 1)
    #pt <- pt[CSmatchpos, ]
      } else { 
    ion.mode <- "Neg"
    meta$baseline <- baselines$neg.bl
    pt <- read.delim("data/EPIC Cross sectional RP NEG Feature table.txt", skip=4, row.names = 1)
    #pt <- pt[CSmatchneg, ]
  }

  #use sampvec to get 451 subjects only
  mat <- t(pt[ , sampvec])
  print(paste(c("Observations:", "Features:"), dim(mat)))

  #replace matrix zeros or NAs with 1s
  mat <- ifelse(mat == 0 | is.na(mat), 1, mat)

  #filter features present in less than 300 samples
  filt <- colSums(mat > 1) > min.sample
  mat  <- mat[ , filt]

  #filter the original peak table (for extraction of feature data)
  pt <- pt[filt, sampvec]

  ### Read and process metadata ------------------------------------------------------------

  #subset samples only
  labs <- meta[meta$present.pos == T & meta$stype == "SA", ]

  #Create food intake object, get quartiles and define intake categories
  foodcol <- labs[, food]
  qfood   <- quantile(foodcol, na.rm = T)
  cats    <- cut(foodcol, qfood, include.lowest = T)

  qmeds <- aggregate(mat, list(classes=cats), median)
  #get features which increase by cof cons quartile (with conditions)
  increasing <- 
      function(x) { which.max(x) == 4 &                         # highest must be quartile 4
                   (which.min(x) == 1 | which.min(x) == 2) &    # lowest must be Q1 or Q2
                    x[3] > x[1] }                               # Q3 must be greater than Q1

  ind <- apply(as.matrix(qmeds[, -1]), 2, increasing)

  #subset increasing only if necessary
  filtmat <- if(incr == T) mat[, ind] else mat

  #impute missings with half minimum value
  library(zoo)
  if(impute == T) filtmat <- na.aggregate(filtmat, FUN = function(x) median(x)/2 )
  
  #subset increasing only and log transform
  logmat <- log(filtmat)
  print(paste("Testing", ncol(logmat), "features..."))
  
  ### Correlation and partial correlation ###-------------------------

  #correlation test food intake with intensities of selected features
  lcor <- apply(logmat, 2, function(x) cor.test(foodcol[x > 0], x[x > 0])  )

  #get raw correlation coefficients
  rcor <- unlist(sapply(lcor, "[", 4))
  rawp <- p.adjust(unlist(sapply(lcor, "[", 3)), method = "BH")
  
  #with purrr
  # rawcor <- map_dfr(lcor, tidy, .id = "feature") %>% mutate(p.valueBH = p.adjust(p.value, method = "BH"))

  #partial correlation controlling for covariates
  #define function to apply to food and intensity data. Change food and covariates as required

  partialcor <- function(x) {
  
    #association with coffee intake. 
    #All obs are passed but logcof automatically has NAs removed and is also subset for nonzero cons
    ## cups
    #xres <- residuals(lm(cups   ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ]))
    ## volume
    mod  <- lm(logcof   ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ])
    xres <- residuals(mod)
    yres <- residuals(lm(x[x > 0] ~ centre + sex + R_BMI + smoke + type.plate + baseline, data = labs[x > 0, ]))
  
    #association with alcohol intake
    #res <- residuals(lm(logQe_Alc ~ country + sex + R_BMI + logcof + smoke, data = labs[x > 0, ]))
    #res <- residuals(lm(x[x > 0]  ~ country + sex + R_BMI + logcof + smoke + type.plate + baseline, data = labs[x > 0, ]))
    #print(mod$call)
    
    #print(paste("Subset size = ", length(xres)))
  
    cor.test(xres, yres)
  }

  #apply function across matrix and extract p-values
  lpcor <- apply(logmat, 2, partialcor)
  Pcor  <- unlist(sapply(lpcor, "[", 4))
  pval  <- p.adjust(unlist(sapply(lpcor, "[", 3)), method = "BH")
  
  # pcor <- map_dfr(lpcor, tidy, .id = "feature") %>% mutate(p.valueBH = p.adjust(p.value, method = "BH"))

  ### Extract feature data ---------------------------------------------------------------------

  #get mass and RT from peak table and add a feature number for joining. First round cols
  #split mass and RT
  splitcl <- rownames_to_column(pt, var = "ID") %>% 
    select(ID) %>%
    separate(ID, into = c("Mass", "RT"), sep = "@", convert = T) %>% 
    mutate(mode = ion.mode, feat = 1:n()) %>% 
    select(mode, Mass, RT, feat) %>%
    mutate(feature = feat)
  
  #Give pos and neg mode data unique feature numbers
  disc <- if(pos == T) splitcl else splitcl %>% mutate(feature = feat + 10000)

  mass.RT <- disc[ind, ]
  med.int <- round(apply(mat[, mass.RT$feat], 2, median))

  #count detections for each feature
  detect  <- apply(mat[, ind], 2, function(x) sum(x > 1))
  #put all the info together, filter by p-value and sort
  disctbl <- data.frame(mass.RT, detect, med.int, rcor, rawp, Pcor, pval) %>%
    filter(pval < pcutoff) %>% arrange(desc(Pcor))

  return(disctbl)
}

### Coffee and alcohol analyses --------------------------------------------------------------------

#1 Coffee, filter increasing, filter not present in three quarters of all samples (n=340)
cofpos <- intake.corr("logcof", incr = T, pcutoff = 0.05)
cofneg <- intake.corr("logcof", incr = T, pos = F, pcutoff = 0.05)

# Coffee, cups, not increasing (for Manhattan set incr to F and pcutoff to 1)
cofpos <- intake.corr("cups", incr = F, impute = F, pcutoff = 1)
cofneg <- intake.corr("cups", incr = F, pos = F, pcutoff = 1)
cof    <- cofpos %>% bind_rows(cofneg) %>% arrange(-Pcor)
#saveRDS(disc.cof, file="Coffee features Manhattan.rds")

#2. Alcohol grams consumed, no increasing filter, missings not imputed:
alcpos <- intake.corr("logQe_Alc", pos = T, incr = F, impute = F, min.sample = 225, pcutoff = 0.05)
alcneg <- intake.corr("logQe_Alc", pos = F, incr = F, impute = F, min.sample = 225, pcutoff = 0.05)
alc   <- alcpos %>% bind_rows(alcneg) %>% arrange(-Pcor)

### Correlation heatmap ----------------------------------------------------------------------

### Read in feature tables ###
ptpos <- read.delim("EPIC Cross sectional RP POS Feature table.txt", skip=4, row.names = 1)
ptneg <- read.delim("EPIC Cross sectional RP NEG Feature table.txt", skip=4, row.names = 1)

ptpos <- ptpos[, sampvec ]
ptneg <- ptneg[, sampvec ]

#merge pos and neg data
disctbl <- rbind(cofpos, cofneg)
saveRDS(disctble, file="coffee_corr_features.rds")

posfilt <- rowSums(ptpos > 1) > 340 #gives logical vector 
negfilt <- rowSums(ptneg > 1) > 340
posmat <- ptpos[ posfilt, ]
negmat <- ptneg[ negfilt, ]

#subset by discriminant feature number
discmatpos <- posmat[ cofpos1$feat,  ] %>% t
discmatneg <- negmat[ cofneg1$feat,  ] %>% t

discmat <- cbind(discmatpos, discmatneg)

#filter the original peak table (for extraction of feature data)
ptpos <- ptpos[cofpos1$feat, ]
ptneg <- ptneg[cofneg1$feat, ]

#log2 transform data and make colnames
logmat   <- log2(discmat)
colnames(logmat) <- paste(disctbl$mode, disctbl$med.int, disctbl$Mass, "@", round(disctbl$RT, 3))
colnames(logmat) <- paste(disctbl$feature)
#colnames(logmat) <- disctbl$feature
cormat   <- cor(logmat)
colnames(cormat) <- NULL

require(corrplot)
#assign to object for cluster order
cplot <- corrplot(cormat, method = "square", tl.col = "black", order = "hclust", tl.cex = 0.5, 
      cl.ratio = 0.2, cl.align = "l", cl.pos = "r", cl.cex = 0.6, mar = c(1,1,1,1))

#Get cluster order to paste into Excel
writeClipboard(rownames(cplot))

### Manhattan plot for coffee updated for new pos and neg data ---------------------------------------
#Load features with partial correlations
library(tidyverse)
disc.cof <- readRDS("Coffee features Manhattan.rds")

#Make colour vector for stripes, 5934 features
colvec   <- c(rep(c(rep("col1", 500), rep("col2", 500)), 5), rep("col1", 500), rep("col2", 434))

#Make data frame for Manhattan ggplot. Added variables: order, the factor of stripe colours, 
#conditional mutate of this for the final colour vector
df <- disc.cof %>% arrange(RT) %>% 
          mutate(order     = 1:5934, 
                 pointcol  = as.factor(colvec), 
                 direction = ifelse(Pcor > 0, "pos", "neg"), 
                 absPcor   = abs(Pcor),
                 signif    = ifelse(pval < 0.05 & direction == "pos", "col3", pointcol),
                 signif2   = ifelse(pval < 0.05 & direction == "neg", "col4", signif),
                 signif3   = ifelse(pval < 0.05, "empty", "filled"),
                 mz        = ifelse(absPcor > 0.30, paste("m/z", Mass, sep = " "), NA),
                 mz2       = ifelse(pval > 0.9953192, paste("m/z", Mass, sep = " "), NA))

#Calculate raw and adjusted cut points for plot
#signif <- filter(df, pval < 0.05)

#Call to ggplot (customise as needed)
#For publication, reduce point size and font size
library(ggplot2)
#ggplot(df, aes(x = order, y = absPcor, colour = signif2)) + 
ggplot(df, aes(x = order, y = -log10(rawp), colour = signif2, shape = signif3)) +
  geom_hline(yintercept = c(0.2, 0.4), linetype = c("solid"), colour = "grey") +
  #geom_hline(yintercept = c(-log10(0.05), -log10(0.9953192)), linetype = "dashed", colour = "black") +
  geom_point(size = 1) + theme_bw(base_size = 10) +
  #scale_colour_manual(values = c("darkgrey", "black", "red", "blue")) +
  scale_colour_manual(values = c("darkgrey", "black", "black", "black")) +
  scale_shape_manual(values = c(1, 16)) + 
  scale_x_continuous(name = "Elution order (increasing lipophilicity)", expand = c(0.02,0)) +
  scale_y_continuous(name = "Partial Pearson correlation coefficient", expand = c(0.01, 0.01)) +
  #scale_y_continuous(name = "-log10(pvalue)", expand = c(0.01, 0.01)) +
  geom_text(aes(label = mz2), hjust = -0.1, vjust = 0, size = 1.5, colour = "black") +
  #geom_hline(yintercept = -0.00, colour = "white", size = 2) +
  theme(legend.position = "none", 
        text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#Plot for publication one-column
ggsave("manhattan coffee paper1.png", width = 100, height = 60, units = "mm")
ggsave("manhattan coffee paper1.svg", width = 100, height = 60, units = "mm")


