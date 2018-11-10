# Match features from metabolomics datasets

# Get common features between CS and HCC datasets and use this as a starting point for analysis

feature.match <- function(RTtol = 0.1, study = c("cs", "hcc")) {

  library(haven)
  library(tidyverse)
  library(fuzzyjoin)
  # read in CS data and HCC data
  
  # CS: 
  cspos <- read_tsv("data/EPIC Cross sectional RP POS Feature table.txt", skip=4) %>% select(1) %>% mutate(CS.feat = 1:n())
  csneg <- read_tsv("data/EPIC Cross sectional RP NEG Feature table.txt", skip=4) %>% select(1) %>% mutate(CS.feat = 1:n())
  
  # HCC
  hccpos <- read_csv("data/EPIC liver cancer 2016 RP POS feature table.csv") %>% 
    select(1) %>% slice(-(83:84)) %>% mutate(HCC.feat = 1:n())
  hccneg <- read_csv("data/EPIC liver cancer 2016 RP Neg Feature Table.csv") %>% 
    select(1) %>% mutate(HCC.feat = 1:n())
  
  # Join pos mode features
  csfeat <- cspos  %>% separate(Compound, into = c("Mass", "RT"), sep = "@", convert = T)
  hccfeat <- hccpos %>% separate(Compound, into = c("Mass", "rt"), sep = "@", convert = T)
  
  posjoin <- difference_inner_join(csfeat, hccfeat, max_dist = 0.005, distance_col = "massdiff") %>% 
    filter(abs(RT - rt) < 0.1) #%>% arrange(CS.feat)
  
  # Join neg mode features
  
  csfeat1 <- csneg  %>% separate(Compound, into = c("Mass", "RT"), sep = "@", convert = T)
  hccfeat1 <- hccneg %>% separate(Compound, into = c("Mass", "rt"), sep = "@", convert = T)
  
  negjoin <- difference_inner_join(csfeat1, hccfeat1, max_dist = 0.005, distance_col = "massdiff") %>% 
    filter(abs(RT - rt) < RTtol) #%>% arrange(CS.feat)
  
  # 2922 features matched in pos and 1243 in neg mode
  
  # Get unique vector of CSS features that matched. These are then used in Intake_correlation to filter starting FTs
  
  if(study == "cs") {
  
    v1 <- unique(posjoin$CS.feat)
    v2 <- unique(negjoin$CS.feat)
  
    } else if (study == "hcc") {
    
    v1 <- unique(posjoin$HCC.feat)
    v2 <- unique(negjoin$HCC.feat)
  
    }
  output <- list(v1, v2)
}  


output1 <- feature.match(study = "cs")
output2 <- feature.match(study = "hcc")

# ----

CSmatchpos <- output1[[1]]
CSmatchneg <- output1[[2]]

HCCmatchpos <- sort(output2[[1]])
HCCmatchneg <- sort(output2[[2]])


# loop to get tolerance vs feature numbers (doesn't work yet)
tol <- seq(0.05, 0.15, by  = 0.005)
output <- numeric(length(tol))
for(i in 1:length(tol)) {
  output[i] <-
  difference_inner_join(csfeat, hccfeat, max_dist = 0.005, distance_col = "massdiff") %>% 
    filter(abs(RT - rt) < tol[i]) %>% nrow()
}

plot(tol, output, type = "o", xlab = "Retention time tolerance", ylab = "No. features matched")
