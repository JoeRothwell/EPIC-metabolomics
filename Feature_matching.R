# Match features from metabolomics datasets

# Get common features between CS and HCC datasets and use this as a starting point for analysis

feature.match <- function(RTtol = 0.1, study = c("cs", "hcc"), mode = c("pos", "neg")) {

  library(haven)
  library(tidyverse)
  library(fuzzyjoin)
  # Get CS data and HCC data (feature names only)
  # Pos or neg mode for CS and HCC
  
  if(mode == "pos") {
  cs <- read_tsv("data/EPIC Cross sectional RP POS Feature table.txt", skip=4) %>% select(1) %>% mutate(CS.feat = 1:n())
  hcc <- read_csv("data/EPIC liver cancer 2016 RP POS feature table.csv") %>% 
    select(1) %>% slice(-(83:84)) %>% mutate(HCC.feat = 1:n())
  
    } else if (mode == "neg") {
      
  cs <- read_tsv("data/EPIC Cross sectional RP NEG Feature table.txt", skip=4) %>% select(1) %>% mutate(CS.feat = 1:n())
  hcc <- read_csv("data/EPIC liver cancer 2016 RP Neg Feature Table.csv") %>% 
    select(1) %>% mutate(HCC.feat = 1:n())
  }
  
  # Join features
  csfeat  <- cs %>% separate(Compound, into = c("Mass", "RT"), sep = "@", convert = T)
  hccfeat <- hcc %>% separate(Compound, into = c("Mass", "rt"), sep = "@", convert = T)
  
  joindf <- difference_inner_join(csfeat, hccfeat, max_dist = 0.005, distance_col = "massdiff") %>% 
    filter(abs(RT - rt) < RTtol) #%>% arrange(CS.feat)
  
  return(joindf)
  
  # 2922 features matched in pos and 1243 in neg mode
  
  # Get unique vector of CSS features that matched. These are then used in Intake_correlation to filter starting FTs
  if(study == "cs") v1 <- unique(joindf$CS.feat) else if (study == "hcc") v1 <- unique(joindf$HCC.feat)
} 

pos.match <- feature.match(mode = "pos")
neg.match <- feature.match(mode = "neg")

# Vectors of matched features
CSmatchpos <- unique(pos.match$CS.feat)
CSmatchneg <- unique(neg.match$CS.feat)
HCCmatchpos <- unique(pos.match$HCC.feat)
HCCmatchneg <- unique(neg.match$HCC.feat)

# ----

# Use vectors to subset peak tables and run partial correlation function
# Now functions from Intake_correlation scripts. Output is assigned alcpos.cs etc

# ----

# Pos mode

alcpos.cs1 <- alcpos.cs %>% select(Pcor : CS.feat)
alcpos.hcc1 <- alcpos.hcc %>% select(Pcor : HCC.feat)
pos.match1 <- pos.match %>% left_join(alcpos.cs1, by = "CS.feat") %>% 
  left_join(alcpos.hcc1, by = "HCC.feat")

# Plot correlations
plot(pos.match1$Pcor.x, pos.match1$Pcor.y, xlab = "Correlation CS", ylab = "Correlation HCC")

# Neg mode

alcneg.cs1  <- alcneg.cs %>% select(Pcor : CS.feat)
alcneg.hcc1 <- alcneg.hcc %>% select(Pcor : HCC.feat)
neg.match1  <- neg.match %>% left_join(alcneg.cs1, by = "CS.feat") %>% 
  left_join(alcneg.hcc1, by = "HCC.feat")

# Plot correlations
plot(neg.match1$Pcor.x, neg.match1$Pcor.y, xlab = "Correlation CS", ylab = "Correlation HCC")

# Plot all 

library(ggplot2)
library(ggrepel)

all.match <- bind_rows(pos.match1, neg.match1, .id = "mode") %>% 
  mutate(mz = ifelse(Pcor.x > 0.2 & Pcor.y > 0.2, paste(Mass.x, round(RT, 2), sep = "@"), NA),
         pointcol= ifelse(Pcor.x > 0.2 & Pcor.y > 0.2, "A", "B"))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

all.match %>%
  ggplot(aes(x = Pcor.x, y = Pcor.y, shape = mode, colour = pointcol)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("black", "darkgrey")) +
  xlab("Partial correlation cross-sectional") + ylab("Partial correlation HCC") +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = mz), size = 3, direction = "both", nudge_x = 0, point.padding = 0.4,
                  segment.size = 0.25, force = 0.2) +
  theme(legend.position = "none") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.5))








# ----

# loop to get tolerance vs feature numbers (doesn't work yet)
tol <- seq(0.05, 0.15, by  = 0.005)
output <- numeric(length(tol))
for(i in 1:length(tol)) {
  output[i] <-
  difference_inner_join(csfeat, hccfeat, max_dist = 0.005, distance_col = "massdiff") %>% 
    filter(abs(RT - rt) < tol[i]) %>% nrow()
}

plot(tol, output, type = "o", xlab = "Retention time tolerance", ylab = "No. features matched")
