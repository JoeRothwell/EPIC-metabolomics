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

# Cross sectional study feature tables
ptpos <- read.delim("EPIC Cross sectional RP POS Feature table.txt", skip=4)
ptneg <- read.delim("EPIC Cross sectional RP NEG Feature table.txt", skip=4)
# HCC
library(tidyverse)
hcpos <- read_csv("EPIC liver cancer 2016 RP POS feature table.csv") %>% slice(-(83:84))
hcneg <- read_csv("EPIC liver cancer 2016 RP NEG feature table.csv")


feature_match <- function(csdat, hccdat, RTtol = 0.1, study = c("cs", "hcc"), filt = NULL) {

  # Get CS data and HCC data (feature names only)
  # Pos or neg mode for CS and HCC
  library(tidyverse)
  cs <- csdat %>% select(1) %>% mutate(CS_feat = 1:n())
  hcc <- hccdat %>% select(1) %>% mutate(HCC_feat = 1:n())
    
  # Filter by the ordered features
  if(!is.null(filt)) cs <- cs[filt, ]
  
  # Join features
  csfeat  <- cs %>% separate(Compound, into = c("Mass", "RT"), sep = "@", convert = T)
  hccfeat <- hcc %>% separate(Compound, into = c("Mass", "rt"), sep = "@", convert = T)

  library(fuzzyjoin)  
  output <- difference_inner_join(csfeat, hccfeat, max_dist = 0.005, distance_col = "massdiff") %>% 
    filter(abs(RT - rt) < RTtol) #%>% arrange(CS.feat)
  
  return(output)
  # 2922 features matched in pos and 1243 in neg mode
  
  # Get unique vector of CSS features that matched. These are then used in Intake_correlation to filter starting FTs
  if(study == "cs") v1 <- unique(joindf$CS_feat) else if (study == "hcc") v1 <- unique(joindf$HCC_feat)
} 

matched_pos1 <- feature_match(ptpos, hcpos)
matched_neg1 <- feature_match(ptpos, hcpos)

# Method 2 Match significant CS features with all HCC

pos_cs <- intakecor_cs(ptpos, food = "Qe_Alc", incr = F, minsamp = 250, pcutoff = 0.05, matchvec = NULL)
neg_cs <- intakecor_cs(ptneg, food = "Qe_Alc", incr = F, minsamp = 250, pcutoff = 0.05, matchvec = NULL)
filt1 <- pos_cs$feat
filt2 <- neg_cs$feat

matched_pos_filt <- CS_HCC_match(mode = "pos", filt = filt1)
matched_neg_filt <- CS_HCC_match(mode = "neg", filt = filt2)