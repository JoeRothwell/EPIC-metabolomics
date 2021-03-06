# Finds features associated with food intake by partial correlation
# Coffee QGE130301; Red meat QGE0701; Fish QGE0801; Fruits, nuts and seeds QGE04; Alc beverages QGE14;
# Beer 140301; Total dietary fibre QEFIBT

# Cross sectional study feature tables
ptpos <- read.delim("EPIC Cross sectional RP POS Feature table.txt", skip=4, row.names = 1)
ptneg <- read.delim("EPIC Cross sectional RP NEG Feature table.txt", skip=4, row.names = 1)
# HCC
hcpos <- read_csv("EPIC liver cancer 2016 RP POS feature table.csv") %>% slice(-(83:84))
hcneg <- read_csv("EPIC liver cancer 2016 RP NEG feature table.csv")

# Function for baseline correction (note: following functions use saved baselines object)
baselines <- function() {
  
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
#bl <- baselines()

# Functions for intake correlation
intakecor_cs <- function(dat, food = "cof", incr = T, impute = F, pcutoff = 0.05, 
                         minsamp = 340, model = T, matchvec = NULL) {
  
  require(tidyverse)
  #see baseline correction.R for details of baseline co-variate
  bl <- readRDS("baselines_pos_neg.rds")
  
  ### Define and process sample metadata (food, lifestyle, technical). Get alcohol (g) and BMI data
  alc_g <- read.csv("alcohol.csv") %>% select(Idepic, Country=country, Qe_Alc, R_BMI)
  cupvol <- data.frame(country = levels(alc_g$Country), cupvol = c(146.59, 209.32, 135.48, 55.2))
  
  # Read in main metadata, add alcohol and BMI data and log transformed intakes
  meta  <- read.csv("cs_metadata.csv") %>% left_join(alc_g, by="Idepic") %>% 
    left_join(cupvol, by = "country") %>% mutate(cups = cof/cupvol)
  
  # Filter the 498 obs to get the 451 measured subjects only (excluding 204 and 355 not properly injected in pos)
  sampvec <- meta$present.pos == T & meta$stype == "SA"
  
  ### Read and process metabolomics data ------------------------------------------------------------
  
  if(nrow(dat) == 5658) ionmode <- "Pos" else ionmode <- "Neg"
  if(nrow(dat) == 5658) meta$baseline <- bl$pos.bl else meta$baseline <- bl$neg.bl
  
  # Add a feature number to the rownames
  rownames(dat) <- paste(rownames(dat), 1:nrow(dat), sep = "@")
  
  # Subset by a vector of features if needed
  if(!is.null(matchvec)) dat <- dat[matchvec, ] else dat
  
  #use sampvec to get 451 subjects only
  mat <- t(dat[ , sampvec])
  print(paste(c("Observations:", "Features:"), dim(mat)))
  
  #replace matrix zeros or NAs with 1s
  mat <- ifelse(mat == 0 | is.na(mat), 1, mat)
  
  # Feature filters ----
  
  # 1. Keep features present in more than n samples only
  filt <- colSums(mat > 1) > minsamp
  mat  <- mat[ , filt]

  # subset 451 samples only (stype = SA)
  labs <- meta %>% filter(present.pos == T & stype == "SA")
  
  # 2. Keep features that are increasing. Get quartiles and define intake categories
  foodcol <- labs[, food]
  cats <- cut_number(foodcol, 4)
  qmeds <- data.frame(cats, mat) %>% group_by(cats) %>% summarise_all(median)
  
  # get features which increase by cof cons quartile (with conditions)
  increasing <- 
    function(x) { which.max(x) == 4 &               # highest must be quartile 4
        (which.min(x) == 1 | which.min(x) == 2) &   # lowest must be Q1 or Q2
        x[3] > x[1] }                               # Q3 must be greater than Q1

  # subset increasing only if necessary  
  ind <- map_lgl(qmeds[, -1], increasing)
  filtmat <- if(incr == T) mat[, ind] else filtmat <- mat
  
  # Get feature metadata
  md <- tibble(medint = apply(filtmat, 2, median), 
               detect = apply(filtmat, 2, function(x) sum(x > 1)))
  
  ### Correlation and partial correlation. First impute and log data
  library(zoo)
  if(impute == T) filtmat <- na.aggregate(filtmat, FUN = function(x) min(x)/2 )
  logmat <- log2(filtmat)
  if(model == F) return(logmat)
  
  # correlation test food intake with intensities of selected features
  print(paste("Testing", ncol(logmat), "features..."))
  print(paste("Correlation with", food))
  lcor <- apply(logmat, 2, function(x) cor.test(foodcol[x > 0], x[x > 0])  )
  
  # partial correlation controlling for covariates
  #define function to apply to food and intensity data. Change food and covariates as required
  partialcor <- function(x) {
    #All obs are passed but logcof automatically has NAs removed and is also subset for nonzero cons
    ## cups or volume
    if(food == "cof") {
      
      #mod1 <- lm(cups   ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ])
      mod1 <- lm(log2(cof + 1) ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ])
      mod2 <- lm(x[x > 0] ~ centre + sex + R_BMI + smoke + type.plate + baseline, data = labs[x > 0, ])
      
    }   else if(food == "Qe_Alc") {
      
      mod1 <- lm(log2(Qe_Alc + 1) ~ country + sex + R_BMI + log2(cof+1) + smoke, data = labs[x > 0, ])
      mod2 <- lm(x[x > 0] ~ country + sex + R_BMI + log2(cof+1) + smoke + type.plate + baseline, data = labs[x > 0, ])
      
    }
    cor.test(residuals(mod1), residuals(mod2))
  }
  
  lpcor <- apply(logmat, 2, partialcor)

  # Tidy and extract correlation coefficients
  library(broom)
  rcor <- map_dfr(lcor, tidy, .id = "feat") %>% separate(feat, sep = "@", into = c("Mass", "RT", "feat"))
  pcor <- map_dfr(lpcor, tidy, .id = "feat") %>% 
    mutate(p.valBH = p.adjust(p.value, method = "BH")) %>% select(estimate, p.value, p.valBH)

  # Final partial correlation is called estimate1  
  output <- bind_cols(rcor, pcor, md) %>% bind_cols(data.frame(t(filtmat))) %>% 
    filter(p.valBH < 0.05)
}

# Old HCC (to be deleted; see Intake_correlation_HCC)
intakecor_hcc <- function(dat, minsamp = 65, matchvec = NULL) {
  # Intensity data and metadata need to first be joined. This is easiest to do by putting lab ID in the same
  # order and then binding by col.
  
  if(nrow(dat) == 4321) ionmode <- "Pos" else ionmode <- "Neg"

  # Subset by a vector of features if needed
  if(!is.null(matchvec)) dat <- dat[matchvec, ] else dat
  
  names <- dat %>% select(Compound) %>% pull
  
  # Subset samples and add feature names
  dat <- dat %>% select(contains("LivCan_")) %>% t
  colnames(dat) <- names
  
  # Split sample ID to get codes for joining to metadata. Remove ID cols
  ints <- dat %>% data.frame %>% rownames_to_column() %>% 
    separate(rowname, sep = "_", into=c("id1", "id2"), convert = T) %>% 
    arrange(id2) %>% select(-(1:2))
  
  # Metadata is from Laura's SAS file
  # As for intensities, split ID to get a code for joining
  library(haven)
  meta <- read_sas("newdataset_missings1.sas7bdat") %>% select(-starts_with("Untarg_")) %>%
    separate(Id_Bma, sep = "_", into=c("id1", "id2"), convert = T) %>% arrange(id2)
  
  # Get vector of controls for subsetting
  ints0 <- ints[meta$Cncr_Caco_Live == 0, ]
  
  #Convert variables of interest to factors
  varlist <- c("Country", "Sex", "Smoke_Stat", "Batch_Rppos", "Batch_Rpneg")
  meta <- meta %>% mutate_at(vars(varlist), as.factor)
  
  #Partial correlation between alcohol intake and each feature
  #Get residuals from lm on alcohol intake: BMI, smoke, sex, batch, country, CC status
  
  # Subsetting: controls only (with suffix 0), matched features only
  
  meta0 <- meta[meta$Cncr_Caco_Live == 0, ]
  filt  <- colSums(ints0 > 1) > minsamp
  
  # Feature match filter will be made here
  ints0 <- ints[meta$Cncr_Caco_Live == 0, filt]
  
  logints <- log2(ints0)
  
  # Subset original feature table for extraction of feature data later
  ft.subset <- dat[meta$Cncr_Caco_Live == 0, filt]
  namesfilt <- names[filt]
  
  print(paste("Testing", ncol(logints), "features..."))
  
  #correlation test food intake with intensities of selected features
  lcor <- apply(logints, 2, function(x) cor.test(log2(meta0$Qe_Alc + 1)[x > 0], x[x > 0]) )
  
  partialcor <- function(x) {
    # Get residuals on log transformed alcohol intake and each feature
    mod1 <- lm(log2(Qe_Alc + 1) ~ Country + Sex + Bmi_C + Smoke_Stat, data = meta0[x > 0, ])
    mod2 <- lm(x[x > 0] ~ Country + Sex + Bmi_C + Batch_Rppos + Smoke_Stat, data = meta0[x > 0, ])
    # Correlate residuals               
    cor.test(residuals(mod1), residuals(mod2))
  }
  
  #Applying the partialcor function across the pos and neg matrices columnwise (argument=2)
  pcor <- apply(logints, 2, partialcor)
  
  #get raw and partial correlation coefficients
  rcor <- unlist(sapply(lcor, "[", 4))
  rawp <- p.adjust(unlist(sapply(lcor, "[", 3)), method = "BH")
  Pcor  <- unlist(sapply(pcor, "[", 4))
  pval  <- p.adjust(unlist(sapply(pcor, "[", 3)), method = "BH")
  
  # Tidying of models
  library(broom)
  rcor <- map_dfr(lcor, tidy, .id = "feat") %>% separate(feat, sep = "@", into = c("Mass", "RT", "feat"))
  pcor <- map_dfr(pcor, tidy, .id = "feat") %>% 
    mutate(p.valBH = p.adjust(p.value, method = "BH")) %>% select(estimate, p.value, p.valBH)
  
  # Split names into mass and RT columns  
  massRT <- data.frame(namesfilt) %>%
    separate(namesfilt, into = c("Mass", "RT"), sep = "@", convert = T) %>% 
    mutate(mode = ionmode, feat = (1:ncol(dat))[filt] )
  
  massRT <- if(ionmode == "Neg") massRT %>% mutate(feat1 = feat + 10000) else massRT %>% mutate(feat1 = feat)
  
  # Get median intensities and count detections for each feature (non-log matrix)
  medint <- round(apply(ints0, 2, median))
  detect <- apply(ints0, 2, function(x) sum(x > 1))
  
  # Add feature match variable and subset if necessary
  if(!is.null(matchvec)) HCC_feat <- matchvec[filt] else matchvec
  
  #put everything into a df
  output <- data.frame(massRT, detect, medint, rcor, rawp, Pcor, pval)
  output <- if(!is.null(matchvec)) cbind(output, HCC_feat) else output
  output <- output %>% filter(pval < 1) %>% arrange(desc(Pcor))
}  

# For CS - HCC matching
CS_HCC_match <- function(RTtol = 0.1, study = c("cs", "hcc"), mode = c("pos", "neg"), filt = NULL) {
  
  library(haven)
  library(tidyverse)
  library(fuzzyjoin)
  # Get CS data and HCC data (feature names only)
  # Pos or neg mode for CS and HCC
  
  if(mode == "pos") {
    cs <- read_tsv("EPIC Cross sectional RP POS Feature table.txt", skip=4) %>% select(1) %>% mutate(CS_feat = 1:n())
    hcc <- read_csv("EPIC liver cancer 2016 RP POS feature table.csv") %>% 
      select(1) %>% slice(-(83:84)) %>% mutate(HCC_feat = 1:n())
    
  } else if (mode == "neg") {
    
    cs <- read_tsv("EPIC Cross sectional RP NEG Feature table.txt", skip=4) %>% select(1) %>% mutate(CS_feat = 1:n())
    hcc <- read_csv("EPIC liver cancer 2016 RP Neg Feature Table.csv") %>% 
      select(1) %>% mutate(HCC_feat = 1:n())
  }
  
  # Filter by the ordered features
  if(!is.null(filt)) cs <- cs[filt, ]
  
  # Join features
  csfeat  <- cs %>% separate(Compound, into = c("Mass", "RT"), sep = "@", convert = T)
  hccfeat <- hcc %>% separate(Compound, into = c("Mass", "rt"), sep = "@", convert = T)
  
  output <- difference_inner_join(csfeat, hccfeat, max_dist = 0.005, distance_col = "massdiff") %>% 
    filter(abs(RT - rt) < RTtol) #%>% arrange(CS.feat)
  
  return(output)
  
  # 2922 features matched in pos and 1243 in neg mode
  
  # Get unique vector of CSS features that matched. These are then used in Intake_correlation to filter starting FTs
  if(study == "cs") v1 <- unique(joindf$CS_feat) else if (study == "hcc") v1 <- unique(joindf$HCC_feat)
} 

# Coffee models for paper
cofpos0 <- intakecor_cs(ptpos, incr = T, pcutoff = 0.05, new = F)
cofneg1 <- intakecor_cs(ptneg, incr = T, pcutoff = 0.05)
# or
output <- lapply(list(ptpos, ptneg), intakecor_cs)


# --------------------

# Old analyses (from previous script; will need updating)

#1 Coffee, filter increasing, filter not present in three quarters of all samples (n=340)
cofpos <- intake.corr("logcof", incr = T, pcutoff = 0.05)
cofneg <- intake.corr("logcof", incr = T, pos = F, pcutoff = 0.05)

# Coffee, cups, not increasing (for Manhattan set incr to F and pcutoff to 1)
# Warning: no longer works and not used in analysis
cofpos <- intake.corr("cups", incr = F, impute = F, pcutoff = 1)
cofneg <- intake.corr("cups", incr = F, pos = F, pcutoff = 1)
cof    <- cofpos %>% bind_rows(cofneg) %>% arrange(-Pcor)
#saveRDS(disc.cof, file="Coffee features Manhattan.rds")

#2. Alcohol grams consumed, no increasing filter, missings not imputed:
# Warning: no longer works, replaced by alcohol_study_untarg.R
alcpos <- intake.corr("logQe_Alc", pos = T, incr = F, impute = F, min.sample = 225, pcutoff = 0.05)
alcneg <- intake.corr("logQe_Alc", pos = F, incr = F, impute = F, min.sample = 225, pcutoff = 0.05)
alc   <- alcpos %>% bind_rows(alcneg) %>% arrange(-Pcor)


