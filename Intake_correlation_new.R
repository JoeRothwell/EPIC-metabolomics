# Finds features associated with food intake by partial correlation
# Coffee QGE130301; Red meat QGE0701; Fish QGE0801; Fruits, nuts and seeds QGE04; Alc beverages QGE14;
# Beer 140301; Total dietary fibre QEFIBT

# Function for baseline correction
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

# ----

intakecor <- function(food = "cof", pos = T, incr = T, impute = F, pcutoff = 0.05, min.sample = 340, model = T){
  
  require(tidyverse)
  #see baseline correction.R for details of baseline co-variate
  baselines <- readRDS("prepdata/baselines pos neg.rds")
  
  ### Define and process sample metadata (food, lifestyle, technical). Get alcohol (g) and BMI data
  alc_g <- read.csv("alcohol/alcohol.csv") %>% select(Idepic, Country=country, Qe_Alc, R_BMI)
  cupvol <- data.frame(country = levels(alc_g$Country), cupvol = c(146.59, 209.32, 135.48, 55.2))
  
  # Read in main metadata, add alcohol and BMI data and log transformed intakes
  meta  <- read.csv("data/cs_metadata.csv") %>% left_join(alc_g, by="Idepic") %>% 
    left_join(cupvol, by = "country") %>% mutate(cups = cof/cupvol)
  
  # Filter the 498 obs to get the 451 measured subjects only (excluding 204 and 355 not properly injected in pos)
  sampvec <- meta$present.pos == T & meta$stype == "SA"
  
  ### Read and process metabolomics data ------------------------------------------------------------
  
  if(pos == T) {  
    ionmode <- "Pos"
    meta$baseline <- baselines$pos.bl
    pt <- read.delim("data/EPIC Cross sectional RP POS Feature table.txt", skip=4, row.names = 1)
    #pt <- pt[CSmatchpos, ]
  } else { 
    ionmode <- "Neg"
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
  
  #subset 451 samples only (stype = SA)
  labs <- meta[meta$present.pos == T & meta$stype == "SA", ]
  
  #Create food intake object, get quartiles and define intake categories
  foodcol <- labs[, food]
  qfood   <- quantile(foodcol, na.rm = T)
  cats    <- cut(foodcol, qfood, include.lowest = T)
  
  qmeds <- aggregate(mat, list(classes=cats), median)
  #get features which increase by cof cons quartile (with conditions)
  increasing <- 
    function(x) { which.max(x) == 4 &               # highest must be quartile 4
        (which.min(x) == 1 | which.min(x) == 2) &   # lowest must be Q1 or Q2
        x[3] > x[1] }                               # Q3 must be greater than Q1
  
  if(incr == T) ind <- apply(as.matrix(qmeds[, -1]), 2, increasing) else ind <- rep(T, ncol(mat))
  
  #subset increasing only if necessary
  #filtmat <- if(incr == T) mat[, ind] else mat
  filtmat <- mat[, ind]
  
  #impute missings with half minimum value
  library(zoo)
  if(impute == T) filtmat <- na.aggregate(filtmat, FUN = function(x) min(x)/2 )
  
  #subset increasing only and log transform
  logmat <- log2(filtmat)
  
  if(model == F) return(logmat)
  
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
    
    #All obs are passed but logcof automatically has NAs removed and is also subset for nonzero cons
    ## cups or volume
    if(food == "cof") {

      #mod1 <- lm(cups   ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ])
    mod1 <- lm(log2(cof + 1)  ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ])
    mod2 <- lm(x[x > 0] ~ centre + sex + R_BMI + smoke + type.plate + baseline, data = labs[x > 0, ])
    
      }   else if(food == "Qe_Alc") {

    mod1 <- lm(log2(Qe_Alc + 1) ~ country + sex + R_BMI + log2(cof+1) + smoke, data = labs[x > 0, ])
    mod2 <- lm(x[x > 0] ~ country + sex + R_BMI + log2(cof+1) + smoke + type.plate + baseline, data = labs[x > 0, ])

      }
    
    cor.test(residuals(mod1), residuals(mod2))
  }
  
  #apply function across matrix and extract p-values
  print(paste("Correlation with", food))
  lpcor <- apply(logmat, 2, partialcor)
  Pcor  <- unlist(sapply(lpcor, "[", 4))
  pval  <- p.adjust(unlist(sapply(lpcor, "[", 3)), method = "BH")
  
  # pcor <- map_dfr(lpcor, tidy, .id = "feature") %>% mutate(p.valueBH = p.adjust(p.value, method = "BH"))
  
  ### Extract feature data ---------------------------------------------------------------------
  
  # get mass and RT from peak table and add a feature number for joining. First round cols
  splitcl <- rownames_to_column(pt, var = "ID") %>% 
    select(ID) %>%
    separate(ID, into = c("Mass", "RT"), sep = "@", convert = T) %>% 
    mutate(mode = ionmode, feat = 1:n()) %>% 
    select(mode, Mass, RT, feat) %>%
    mutate(feature = feat)
  
  # Give pos and neg mode data unique feature numbers. Neg numbers start from 10000
  disc <- if(pos == T) splitcl else splitcl %>% mutate(feature = feat + 10000)
  
  massRT <- disc[ind, ]
  
  # Get median intensities and count detections for each feature
  medint <- round(apply(mat[, massRT$feat], 2, median))
  detect  <- apply(mat[, ind], 2, function(x) sum(x > 1))
  
  #put everything into a df
  disctbl <- data.frame(massRT, detect, medint, rcor, rawp, Pcor, pval) %>%
    filter(pval < pcutoff) %>% arrange(desc(Pcor))

}

# Coffee. Set pcutoff = 1 for Manhattan
cofpos <- intakecor(incr = T, pos = T, pcutoff = 0.05)
cofneg <- intakecor(incr = T, pos = F, pcutoff = 0.05)

# Alcohol
alcpos <- intakecor(food = "Qe_Alc", incr = F)
alcneg <- intakecor(food = "Qe_Alc", incr = F, pos = F)

# Filtered by CS/HCC matching
alcpos1 <- intakecor(food = "Qe_Alc", incr = F, min.sample = 250)
alcneg1 <- intakecor(food = "Qe_Alc", incr = F, pos = F, min.sample = 250)

# Coffee matrix only
logmat <- intakecor(incr = F, impute = T, model = F)
scalemat <- scale(logmat)

# Correlation heatmap ----

corrdata <- function(posdisc, negdisc) {
  ptpos <- read.delim("data/EPIC Cross sectional RP POS Feature table.txt", skip=4, row.names = 1)
  ptneg <- read.delim("data/EPIC Cross sectional RP NEG Feature table.txt", skip=4, row.names = 1)
  meta  <- read.csv("data/cs_metadata.csv")
  
  # Subset samples and features by filtering
  samples <- meta$present.pos == T & meta$stype == "SA" 
  posfilt <- rowSums(ptpos > 1) > 340
  negfilt <- rowSums(ptneg > 1) > 340
  ptpos <- ptpos[posfilt, samples ]
  ptneg <- ptneg[negfilt, samples ]
  
  # merge pos and neg data
  disctbl <- rbind(posdisc, negdisc)
  # saveRDS(disctbl, file="coffee_corr_features.rds")
  
  #subset by discriminant feature number
  discmatpos <- ptpos[ posdisc$feat,  ] %>% t
  discmatneg <- ptneg[ negdisc$feat,  ] %>% t
  
  discmat <- cbind(discmatpos, discmatneg)
  
  #filter the original peak table (for extraction of feature data)
  ptpos <- ptpos[posdisc$feat, ]
  ptneg <- ptneg[negdisc$feat, ]
  
  #log2 transform data and make colnames
  logmat   <- log2(discmat)
  colnames(logmat) <- paste(disctbl$mode, disctbl$medint, disctbl$Mass, "@", round(disctbl$RT, 3))
  #colnames(logmat) <- paste(disctbl$feature)
  #colnames(logmat) <- disctbl$feature
  cormat   <- cor(logmat)
  colnames(cormat) <- NULL
  return(cormat)
}
cormat <- corrdata(cofpos, cofneg)

library(corrplot)
cplot <- corrplot(cormat, method = "square", tl.col = "black", order = "hclust", 
        tl.cex = 0.5, cl.ratio = 0.2, cl.align = "l", cl.pos = "r", 
        cl.cex = 0.6, mar = c(1,1,1,1))

#Get cluster order to paste into Excel
writeClipboard(rownames(cplot))

# Manhattan plot ----

manhattandata <- function() {
  library(tidyverse)
  cof <- readRDS("prepdata/Coffee_features_Manhattan.rds")
  
  # Vector of colours for stripes, length 5934
  colvec <- rep(c("col1", "col2"), 6, each = 500)[1:nrow(cof)]
  colvec2 <- rep(brewer.pal(10, "Paired"), each = 600)[1:nrow(cof)]
  
  #Make data frame for Manhattan ggplot. Added variables: order, the factor of colours, 
  #conditional mutate of this for the final colour vector
  df <- cof %>% arrange(RT) %>% 
    mutate(order     = 1:n(), 
           pointcol = colvec2, 
           direction = ifelse(Pcor > 0, "pos", "neg"), 
           signif    = ifelse(pval < 0.05 & direction == "pos", "#000000", pointcol),
           signif2   = ifelse(pval < 0.05 & direction == "neg", "#000000", signif),
           signif3   = ifelse(pval < 0.05, "empty", "filled"),
           #mz2       = ifelse(pval > 0.9953192, paste("m/z", Mass, sep = " "), NA)
           mz        = ifelse(abs(Pcor) > 0.30, paste("m/z", Mass, sep = " "), NA))
}
df <- manhattandata()

#Calculate raw and adjusted cut points for plot
#signif <- filter(df, pval < 0.05)

#For publication, reduce point size and font size
library(ggplot2)
ggplot(df, aes(x = order, y = abs(Pcor), colour = signif2, shape = signif3)) + 
#ggplot(df, aes(x = order, y = -log10(rawp), colour = signif2, shape = signif3)) +
  geom_hline(yintercept = c(0.2, 0.4), linetype = c("dashed"), colour = "grey") +
  #geom_hline(yintercept = c(-log10(0.05), -log10(0.9953192)), linetype = "dashed", colour = "black") +
  geom_point(size = 1) + theme_bw(base_size = 10) +
  #scale_colour_manual(values = c("darkgrey", "black", "red", "blue")) +
  #scale_colour_manual(values = c("darkgrey", "black", "black", "black")) +
  scale_shape_manual(values = c(1, 16)) + 
  scale_x_continuous(name = "Elution order (increasing lipophilicity)", expand = c(0.02,0)) +
  scale_y_continuous(name = "Partial Pearson correlation coefficient", expand = c(0.01, 0.01)) +
  #scale_y_continuous(name = "-log10(pvalue)", expand = c(0.01, 0.01)) +
  geom_text(aes(label = mz), hjust = -0.1, vjust = 0, size = 1.5, colour = "black") +
  #geom_hline(yintercept = -0.00, colour = "white", size = 2) +
  theme(legend.position = "none", 
        text = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Plot for publication one-column
ggsave("manhattan coffee paper1.png", width = 100, height = 60, units = "mm")
ggsave("manhattan coffee paper1.svg", width = 100, height = 60, units = "mm")
