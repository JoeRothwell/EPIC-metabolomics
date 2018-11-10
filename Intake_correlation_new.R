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
    pt <- pt[CSmatchpos, ]
  } else { 
    ionmode <- "Neg"
    meta$baseline <- baselines$neg.bl
    pt <- read.delim("data/EPIC Cross sectional RP NEG Feature table.txt", skip=4, row.names = 1)
    pt <- pt[CSmatchneg, ]
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
  
  #get mass and RT from peak table and add a feature number for joining. First round cols
  splitcl <- rownames_to_column(pt, var = "ID") %>% 
    select(ID) %>%
    separate(ID, into = c("Mass", "RT"), sep = "@", convert = T) %>% 
    mutate(mode = ionmode, feat = 1:n()) %>% 
    select(mode, Mass, RT, feat) %>%
    mutate(feature = feat)
  
  #Give pos and neg mode data unique feature numbers. Neg numbers start from 10000
  disc <- if(pos == T) splitcl else splitcl %>% mutate(feature = feat + 10000)
  
  massRT <- disc[ind, ]
  medint <- round(apply(mat[, massRT$feat], 2, median))
  
  #count detections for each feature
  detect  <- apply(mat[, ind], 2, function(x) sum(x > 1))
  
  #put everything into a df
  disctbl <- data.frame(massRT, detect, medint, rcor, rawp, Pcor, pval) %>%
    filter(pval < pcutoff) %>% arrange(desc(Pcor))

}

# Coffee
cofpos <- intakecor()
cofneg <- intakecor(pos = F)

# Alcohol
alcpos <- intakecor(food = "Qe_Alc", incr = F)
alcneg <- intakecor(food = "Qe_Alc", incr = F, pos = F)

# Filtered by CS/HCC matching
alcpos1 <- intakecor(food = "Qe_Alc", incr = F, min.sample = 250)
alcneg1 <- intakecor(food = "Qe_Alc", incr = F, pos = F, min.sample = 250)

# Coffee matrix only
logmat <- intakecor(incr = F, impute = T, model = F)
scalemat <- scale(logmat)
