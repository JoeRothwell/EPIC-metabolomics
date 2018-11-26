# All functions and workflow for alcohol study

# Intake correlation for CS and HCC study and match function

intakecor_cs <- function(food = "cof", pos = T, incr = T, impute = F, pcutoff = 0.05, min.sample = 340, model = T, matchvec = NULL){
  
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
  
  # Subset by a vector of features if needed
  if(!is.null(matchvec)) pt <- pt[matchvec, ] else pt
  
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
      mod1 <- lm(log2(cof + 1) ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ])
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
  detect <- apply(mat[, ind], 2, function(x) sum(x > 1))
  
  # Get IDs of matched features
  #CS.feat <- if(pos == T) CSmatchpos[filt] else CSmatchneg[filt]
  if(!is.null(matchvec)) CS_feat <- matchvec[filt] else matchvec
  
  #put everything into a df
  output <- data.frame(massRT, detect, medint, rcor, rawp, Pcor, pval)
  output <- if(!is.null(matchvec)) cbind(output, CS_feat) else output
  output <- output %>% filter(pval < pcutoff) %>% arrange(desc(Pcor))
  
}
intakecor_hcc <- function(pos = T, matchvec = NULL) {
  # Intensity data and metadata need to first be joined. This is easiest to do by putting lab ID in the same
  # order and then binding by col.
  
  # Intensity data
  library(tidyverse)
  if(pos == T) {
    ionmode <- "Pos"
    ft <- read_csv("data/EPIC liver cancer 2016 RP POS feature table.csv") %>% slice(-(83:84))
    #ft <- ft[HCCmatchpos, ]
  } else {
    ionmode <- "Neg"
    ft <- read_csv("data/EPIC liver cancer 2016 RP NEG feature table.csv")
    #ft <- ft[HCCmatchneg, ]
  }
  
  # Subset by a vector of features if needed
  if(!is.null(matchvec)) ft <- ft[matchvec, ] else ft
  
  names <- ft %>% select(Compound) %>% pull
  
  # Subset samples and add feature names
  ft <- ft %>% select(contains("LivCan_")) %>% t
  colnames(ft) <- names
  
  # Split sample ID to get codes for joining to metadata. Remove ID cols
  ints <- ft %>% data.frame %>% rownames_to_column() %>% 
    separate(rowname, sep = "_", into=c("id1", "id2"), convert = T) %>% 
    arrange(id2) %>% select(-(1:2))
  
  # Metadata is from Laura's SAS file
  # As for intensities, split ID to get a code for joining
  library(haven)
  meta <- read_sas("data/newdataset_missings1.sas7bdat") %>% select(-starts_with("Untarg_")) %>%
    separate(Id_Bma, sep = "_", into=c("id1", "id2"), convert = T) %>% arrange(id2)
  
  # Get vector of controls for subsetting
  controls <- meta$Cncr_Caco_Live == 0
  ints0 <- ints[controls, ]
  
  #Convert variables of interest to factors
  varlist <- c("Country", "Sex", "Smoke_Stat", "Batch_Rppos", "Batch_Rpneg")
  meta <- meta %>% mutate_at(vars(varlist), as.factor)
  
  #Partial correlation between alcohol intake and each feature
  #Get residuals from lm on alcohol intake: BMI, smoke, sex, batch, country, CC status
  
  # Subsetting: controls only (with suffix 0), matched features only
  
  meta0 <- meta[controls, ]
  filt  <- colSums(ints0 > 1) > 65
  
  # Feature match filter will be made here
  ints0 <- ints[controls, filt]
  
  logints <- log2(ints0)
  
  # Subset original feature table for extraction of feature data later
  ft <- ft[controls, filt]
  namesfilt <- names[filt]
  
  print(paste("Testing", ncol(logints), "features..."))
  
  #correlation test food intake with intensities of selected features
  lcor <- apply(logints, 2, function(x) cor.test(log2(meta0$Qe_Alc + 1)[x > 0], x[x > 0]) )
  
  #get raw correlation coefficients
  rcor <- unlist(sapply(lcor, "[", 4))
  rawp <- p.adjust(unlist(sapply(lcor, "[", 3)), method = "BH")
  
  partialcor <- function(x) {
    
    # Get residuals on log transformed alcohol intake and each feature
    mod1 <- lm(log2(Qe_Alc + 1) ~ Country + Sex + Bmi_C + Smoke_Stat, data = meta0[x > 0, ])
    mod2 <- lm(x[x > 0] ~ Country + Sex + Bmi_C + Batch_Rppos + Smoke_Stat, data = meta0[x > 0, ])
    
    # Correlate residuals               
    cor.test(residuals(mod1), residuals(mod2))
  }
  
  #Applying the partialcor function across the pos and neg matrices columnwise (argument=2)
  lpcorpos <- apply(logints, 2, partialcor)
  
  #extract p-values and correlations (pos mode)
  Pcor  <- unlist(sapply(lpcorpos, "[", 4))
  pval  <- p.adjust(unlist(sapply(lpcorpos, "[", 3)), method = "BH")
  
  massRT <- data.frame(namesfilt) %>%
    separate(namesfilt, into = c("Mass", "RT"), sep = "@", convert = T) %>% 
    mutate(mode = ionmode, feat = 1:n()) %>%
    mutate(feature = feat)
  
  # Get median intensities and count detections for each feature (non-log matrix)
  medint <- round(apply(ints0[, massRT$feat], 2, median))
  detect <- apply(ints0, 2, function(x) sum(x > 1))
  
  # Add feature match variable and subset if necessary
  if(!is.null(matchvec)) HCC_feat <- matchvec[filt] else matchvec
  
  #put everything into a df
  output <- data.frame(massRT, detect, medint, rcor, rawp, Pcor, pval)
  output <- if(!is.null(matchvec)) cbind(output, HCC_feat) else output
  output <- output %>% filter(pval < 1) %>% arrange(desc(Pcor))
}
CS_HCC_match <- function(RTtol = 0.1, study = c("cs", "hcc"), mode = c("pos", "neg"), filt = NULL) {
  
  library(haven)
  library(tidyverse)
  library(fuzzyjoin)
  # Get CS data and HCC data (feature names only)
  # Pos or neg mode for CS and HCC
  
  if(mode == "pos") {
    cs <- read_tsv("data/EPIC Cross sectional RP POS Feature table.txt", skip=4) %>% select(1) %>% mutate(CS_feat = 1:n())
    hcc <- read_csv("data/EPIC liver cancer 2016 RP POS feature table.csv") %>% 
      select(1) %>% slice(-(83:84)) %>% mutate(HCC_feat = 1:n())
    
  } else if (mode == "neg") {
    
    cs <- read_tsv("data/EPIC Cross sectional RP NEG Feature table.txt", skip=4) %>% select(1) %>% mutate(CS_feat = 1:n())
    hcc <- read_csv("data/EPIC liver cancer 2016 RP Neg Feature Table.csv") %>% 
      select(1) %>% mutate(HCC_feat = 1:n())
  }
  
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



# Approach 1: match HCC and CS features and perform separate correlations on the overlap ----

matched_pos <- CS_HCC_match(mode = "pos")
matched_neg <- CS_HCC_match(mode = "neg")

# Extract vectors of matched features
CSpos <- unique(matched_pos$CS_feat)
CSneg <- unique(matched_neg$CS_feat)
HCCpos <- unique(matched_pos$HCC_feat)
HCCneg <- unique(matched_neg$HCC_feat)

# Run correlations on matched features for both CS and HCC
alcpos_cs <- intakecor_cs(food = "Qe_Alc", incr = F, min.sample = 250, pcutoff = 1, matchvec = CSpos)
alcneg_cs <- intakecor_cs(food = "Qe_Alc", incr = F, pos = F, min.sample = 250, pcutoff = 1, matchvec = CSneg)
alcpos_hcc <- intakecor_hcc(matchvec = HCCpos)
alcneg_hcc <- intakecor_hcc(pos = F, matchvec = HCCneg)

# Join pos and neg data together for each study
pos_cs <- alcpos_cs %>% select(Pcor:CS_feat)
pos_hcc <- alcpos_hcc %>% select(Pcor:HCC_feat)
matched_pos <- matched_pos %>% left_join(pos_cs, by = "CS_feat") %>% left_join(pos_hcc, by = "HCC_feat")

neg_cs <- alcneg_cs %>% select(Pcor:CS_feat)
neg_hcc <- alcneg_hcc %>% select(Pcor:HCC_feat)
matched_neg <- matched_neg %>% left_join(neg_cs, by = "CS_feat") %>% left_join(neg_hcc, by = "HCC_feat")

# Plot correlations (4165 in both modes)
plot(matched_pos$Pcor.x, matched_pos$Pcor.y, xlab = "Correlation CS", ylab = "Correlation HCC")
points(matched_neg$Pcor.x, matched_neg$Pcor.y, xlab = "Correlation CS", ylab = "Correlation HCC")
# See Feature_matching for ggplot plot




# Approach 2: Match significant CS features with all HCC ----

# Get significant CS features
alcpos_cs <- intakecor_cs(food = "Qe_Alc", incr = F, min.sample = 250, pcutoff = 0.05, matchvec = NULL)
alcneg_cs <- intakecor_cs(food = "Qe_Alc", incr = F, pos = F, min.sample = 250, pcutoff = 0.05, matchvec = NULL)
filt1 <- alcpos_cs$feat
filt2 <- alcneg_cs$feat

# Match with all HCC features and get vector of matched features
matched_pos_filt <- CS_HCC_match(mode = "pos", filt = filt1)
matched_neg_filt <- CS_HCC_match(mode = "neg", filt = filt2)
f1 <- matched_pos_filt$HCC_feat
f2 <- matched_neg_filt$HCC_feat

pos_hcc <- intakecor_hcc(matchvec = f1)
neg_hcc <- intakecor_hcc(pos = F, matchvec = f2)
