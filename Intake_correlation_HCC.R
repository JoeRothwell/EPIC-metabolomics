# Partial correlation between alcohol intake of HCC controls and spectral features
# New version of Liver CC alcohol model.R 
# matchvec allows the passing of a vector that subsets the initial peak table

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

alcpos.hcc <- intakecor.hcc()
alcneg.hcc <- intakecor.hcc(pos = F)


# Get the subset of HCC features matched with significant alcohol CS features (from Intake_correlation_new.R) 
f1 <- matched_pos_filt$HCC.feat
f2 <- matched_neg_filt$HCC.feat

alcpos.hcc1 <- intakecor_hcc(matchvec = f1)
alcneg.hcc1 <- intakecor_hcc(pos = F, matchvec = f2)



