# Partial correlation between alcohol intake of HCC controls and spectral features
# New version of Liver CC alcohol model.R 

hcpos <- read_csv("EPIC liver cancer 2016 RP POS feature table.csv") %>% slice(-(83:84))
hcneg <- read_csv("EPIC liver cancer 2016 RP NEG feature table.csv")

intakecor_hcc <- function(dat, pos = T, matchvec = NULL) {
  
  # Feature tables
  library(tidyverse)
  if(nrow(dat) == 4321) ionmode <- "Pos" else ionmode <- "Neg"
  names <- dat$Compound
  
  # Subset samples and add feature names
  dat <- dat %>% select(contains("LivCan_")) %>% t
  colnames(dat) <- paste(names, 1:ncol(dat), sep = "@")
  
  # Split sample ID and order by second part
  id <- data.frame(val = rownames(dat)) %>% separate(val, sep = "_", into=c("id1", "id2"), convert = T)
  ints <- cbind(id, dat) %>% arrange(id2) %>% select(-(1:2))
  
  # Metadata from Laura. Split ID to put in same order as ints
  library(haven)
  meta <- read_sas("newdataset_missings1.sas7bdat") %>% select(-starts_with("Untarg_")) %>%
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
  
  # Subset and log transform matrix
  ints0 <- ints[controls, filt]
  
  # Get median intensity and numbers of detections (untransformed)
  md <- tibble(medint = apply(ints0, 2, median), 
               detect = apply(ints0, 2, function(x) sum(x > 1)))
  
  logints <- log2(ints0)
  
  # Subset original feature table for extraction of feature data later
  dat <- dat[controls, filt]
  namesfilt <- names[filt]
  
  print(paste("Testing", ncol(logints), "features..."))
  
  # correlation test food intake with intensities of selected features
  lcor <- apply(logints, 2, function(x) cor.test(log2(meta0$Qe_Alc + 1)[x > 0], x[x > 0]) )
  
  partialcor <- function(x) {
    # Get residuals on log transformed alcohol intake and each feature
    mod1 <- lm(log2(Qe_Alc + 1) ~ Country + Sex + Bmi_C + Smoke_Stat, data = meta0[x > 0, ])
    mod2 <- lm(x[x > 0] ~ Country + Sex + Bmi_C + Batch_Rppos + Smoke_Stat, data = meta0[x > 0, ])
    # Correlate residuals               
    cor.test(residuals(mod1), residuals(mod2))
  }
  
  # Applying the partialcor function across the pos and neg matrices columnwise (argument=2)
  pcor <- apply(logints, 2, partialcor)

  library(broom)
  rcor <- map_dfr(lcor, tidy, .id = "feat") %>% separate(feat, sep = "@", into = c("Mass", "RT", "feat"))
  pcor <- map_dfr(pcor, tidy, .id = "feat") %>% 
    mutate(p.valBH = p.adjust(p.value, method = "BH")) %>% select(estimate, p.value, p.valBH)
  
  # Final partial correlation is called estimate1  
  output <- bind_cols(rcor, pcor, md) %>% bind_cols(data.frame(t(ints0))) %>%
    arrange(desc(estimate1))
}

alcpos.hcc <- intakecor_hcc(hcpos)
alcneg.hcc0 <- intakecor_hcc(hcneg)


# Get the subset of HCC features matched with significant alcohol CS features (from Intake_correlation_new.R) 
f1 <- matched_pos_filt$HCC.feat
f2 <- matched_neg_filt$HCC.feat

alcpos.hcc1 <- intakecor_hcc(matchvec = f1)
alcneg.hcc1 <- intakecor_hcc(pos = F, matchvec = f2)



