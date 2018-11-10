#Performs discriminant analysis by partial correlation
#Coffee QGE130301; Red meat QGE0701; Fish QGE0801; Fruits, nuts and seeds QGE04; Alc beverages QGE14;
#Beer 140301; Total dietary fibre QEFIBT

intake.corr <- function(food, pos = T, incr = T, impute = F, min.sample = 340, pcutoff = 0.05) {
  
  require(tidyverse)
  #see baseline correction.R for details of baseline co-variate
  baselines <- readRDS("prepdata/baselines pos neg.rds")
  
  ### Define and process sample metadata (food, lifestyle, technical) ###
  #Get alcohol (g) and BMI data
  alc_g <- read.csv("alcohol/alcohol.csv") %>% select(Idepic, Qe_Alc, R_BMI)
  
  # Read in main metadata, add alcohol and BMI data and log transformed intakes
  meta  <- read.csv("data/cs_metadata.csv") %>%
           left_join(alc_g, by="Idepic") %>%
           mutate(logcof = log(cof + 1), logalc = log(alcbev + 1), logredmeat = log(redmeat + 1),
           logprocmeat = log(procmeat + 1), logfish = log(fish + 1), logQe_Alc = log(Qe_Alc + 1))
  
  #add cup volumes in ml: France, 146.59; Italy, 55.2; Greece, 135.48; Germany, 209.32
  cupvols <- data.frame(country = levels(meta$country), cupvol = c(146.59, 209.32, 135.48, 55.2))
  meta    <- meta %>% left_join(cupvols, by="country") %>% mutate(cups = cof/cupvol)
  
  
  
  #---------------------------------------------------
  
  ### Define and process metabolomics data ###
  if(pos == T) {     ion.mode <- "Pos"
                     meta$baseline <- baselines$pos.bl 
                     pt <- read.delim("data/EPIC Cross sectional RP POS Feature table.txt", skip=4)
    } else { 
                     ion.mode <- "Neg"
                     meta$baseline <- baselines$neg.bl
                     pt <- read.delim("data/EPIC Cross sectional RP NEG Feature table.txt", skip=4)
    }
  
  #remove compound column and empty data file 186
  mat <- t(pt[ , -c(1, 187)])
  print(paste(c("Observations:", "Features:"), dim(mat)))
  
  
  
  
  
  #replace matrix zeros or NAs with 1s
  mat <- ifelse(mat == 0 | is.na(mat), 1, mat)
  
  #filter features present in less than 300 samples (+ 40 QCs)
  filt <- colSums(mat > 1) > min.sample
  mat  <- mat[ , filt]
  
  #filter the original peak table (for extraction of feature data)
  pt <- pt[filt, ]
  
  #----------------------------------------------------
  
  #removes the one sample missing in pos (gives 497 obs)
  labs <- meta[meta$present == T, ]
  #join fibreffq to allffq and allffq to other metadata
  #allffq   <- read.csv("allffq.csv")
  #fibreffq <- read.csv("Fibers_FFQs.csv")
  #ffq.fib  <- allffq %>% left_join(fibreffq, by="Idepic")
  #labs     <- labs %>% left_join(ffq.fib, by="Idepic")
  
  #Create food intake object, get quartiles and define intake categories
  foodcol <- labs[, food]
  qfood   <- quantile(foodcol, na.rm = T)
  cats    <- cut(foodcol, qfood, include.lowest = T)
  
  qmeds <- aggregate(mat, list(classes=cats), median)
  #get features which increase by cof cons quartile (with conditions)
  increasing <- function(x) { which.max(x) == 4 &     #highest must be quartile 4
      (which.min(x) == 1 | which.min(x) == 2) &       #lowest must be Q1 or Q2
      x[3] > x[1] }                                   #Q3 must be greater than Q1
  
  if(incr == T) ind <- apply(as.matrix(qmeds[, -1]), 2, increasing) else ind <- rep(T, ncol(mat))
  
  #subset increasing only if necessary
  filtmat <- mat[, ind]
  
  #impute missings with half minimum value
  if(impute == T) {
    filtmat <- apply(filtmat, 2, function(x) { x[which(x == 1)] <- min(x[x > 1], na.rm = T) / 2
    return(x) } )
  } else filtmat 
  
  #subset increasing only and log transform
  logmat <- log(filtmat)
  print(paste("Testing", ncol(logmat), "features..."))
  
  #--------------------------------------------------------------
  
  ### Correlation and partial correlation ###
  
  #correlation test food intake with intensities of selected features
  lcor <- apply(logmat, 2, function(x) cor.test(foodcol[x > 0], x[x > 0])  )
  
  #get raw correlation coefficients
  rcor <- unlist(sapply(lcor, "[", 4))
  rawp <- p.adjust(unlist(sapply(lcor, "[", 3)), method = "BH")
  
  #partial correlation controlling for covariates
  #define function to apply to food and intensity data. Change food and covariates as required
  
  
  

  partialcor <- function(x) {
    
    #association with coffee intake. 
    #All obs are passed but logcof automatically has NAs removed and is also subset for nonzero cons
    ## cups
    #xres <- residuals(lm(cups   ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ]))
    ## volume
    #xres <- residuals(lm(logcof   ~ centre + sex + R_BMI + smoke, data = labs[x > 0, ]))
    #yres <- residuals(lm(x[x > 0] ~ centre + sex + R_BMI + smoke + type.plate + baseline, data = labs[x > 0, ]))
    
    
    #association with alcohol intake
    xres <- residuals(lm(logQe_Alc ~ country + sex + R_BMI + logcof + smoke, data = labs[x > 0, ]))
    yres <- residuals(lm(x[x > 0]  ~ country + sex + R_BMI + logcof + smoke + type.plate + baseline, data = labs[x > 0, ]))

    
    print(paste("Subset size = ", length(xres)))
    
    cor.test(xres, yres)
  }
  
  #apply function across matrix and extract p-values
  lpcor <- apply(logmat, 2, partialcor)
  Pcor  <- unlist(sapply(lpcor, "[", 4))
  pval  <- p.adjust(unlist(sapply(lpcor, "[", 3)), method = "BH")
  
  #-----------------------------------------------------------------
  
  ### Extract feature data ###
  
  #get mass and RT from peak table and add a feature number for joining. First round cols
  #split mass and RT
  splitcl <- pt %>% separate(Compound, into = c("Mass", "RT"), sep = "@", convert = T) %>% 
      mutate(mode = ion.mode, feature = c(1:nrow(pt))) %>% 
      select(mode, feature, Mass, RT)
  
  mass.RT <- splitcl[ind, ]
  med.int <- round(apply(mat[, mass.RT$feature], 2, median))
  
  
  
  
  
  #count detections for each feature
  #get vector of samples only without QCs and blanks
  samples <- !is.na(labs$country)
  detect  <- apply(mat[samples, ind], 2, function(x) sum(x > 1))
  #put all the info together, filter by p-value and sort
  disctbl <- data.frame(mass.RT, detect, med.int, rcor, rawp, Pcor, pval) %>%
    filter(pval < pcutoff) %>% arrange(desc(Pcor))
  
  return(disctbl)
}

#Analyses for different projects

#1. Coffee, filter increasing, filter not present in two thirds of samples
disctbl.pos1 <- intake.corr("logcof", incr = T, impute = F, pcutoff = 0.05)
disctbl.neg1 <- intake.corr("logcof", incr = T, pos = F, pcutoff = 0.05)

# Coffee, cups, not increasing
disctbl.pos <- intake.corr("cups", incr = T, impute = F, pcutoff = 0.05)
disctbl.neg <- intake.corr("cups", incr = T, pos = F, pcutoff = 0.05)

#2. Alcohol grams consumed, no increasing filter, missings not imputed: (for Manhattan set pcutoff = 1)
posdisc <- intake.corr("logQe_Alc", pos = T, incr = F, impute = F, min.sample = 265, pcutoff = 0.05)
negdisc <- intake.corr("logQe_Alc", pos = F, incr = F, impute = F, min.sample = 265, pcutoff = 0.05)
disc <- posdisc %>% bind_rows(negdisc) %>% arrange(-Pcor)
#saveRDS(alc.all.ss, "Alcohol features Manhattan.rds")

#3. Alcohol grams consumed, no increasing filter, missings imputed:
disc.alc.pos2 <- intake.corr("logQe_Alc", pos = T, incr = F, impute = T, min.sample = 265)
disc.alc.neg2 <- intake.corr("logQe_Alc", pos = F, incr = F, impute = T, min.sample = 265)
alc.all.imp <- disc.alc.pos2 %>% bind_rows(disc.alc.neg2) %>% arrange(-Pcor)

#bind pos and neg data together and write to csv
allcor <- bind_rows(disctbl.pos, disctbl.neg) %>% arrange(desc(Pcor))
write.csv(allcor, "Correlated coffee intake cups.csv")

#-------------------------------------------------------------

### Correlation heatmap ### 
#First retrieve intensities, log transform and generate cor matrix

discheatmap <- function(disctbl, pos = T, min.sample = 340) {
  
  require(tidyverse)
  ### Define and process metabolomics data ###
  if(pos == T) { ptpos <- read.csv("EPIC RP POS Corrected final peak table into MPP.csv")
                 mat   <- t(ptpos[, -c(1:5)])
                 ptpos$RT <- round(ptpos$RT, 3)
    } else { 
                 ptneg <- read.delim("EPIC serum RP NEG corrected.txt", header=T, skip=4)  
                 mat   <- t(ptneg[, -1]) 
    }
  
  mat      <- ifelse(mat == 0 | is.na(mat), 1, mat)
  
  #filter features present in less than 300 samples (+ 40 QCs)
  filt <- colSums(mat > 1) > min.sample
  mat <- mat[ , filt]
  
  #filter the original peak table (for extraction of feature data)
  if(pos == T) ptpos <- ptpos[filt, ] else ptneg <- ptneg[filt, ]
  
  discmat  <- mat[ , disctbl$feature]
  logmat   <- log2(mat)[, disctbl$feature]
  cormat   <- cor(logmat)
  
  if(pos == T) {
    rownames(cormat) <- paste(round(apply(discmat, 2, median)),",", ptpos[disctbl$feature, 1],"@",
                              ptpos[disctbl$feature, 2])
    } else {
    rownames(cormat) <- ptneg[disctbl$feature, 1]
    }
  
  require(corrplot)
  corrplot(cormat, method="square", tl.col="black", order="hclust", tl.cex=0.75, 
           cl.ratio=0.2, cl.align="l", cl.pos="r", cl.cex=0.6, mar=c(1,1,1,1))
}

discheatmap(disctbl)
discheatmap(disctbl.neg, pos = F)

#updated for new pos and neg data
  
  require(tidyverse)
  ### Define and process metabolomics data ###
  ptpos <- read.delim("EPIC Cross sectional RP POS Feature table.txt", skip=4, row.names = 1)
  ptneg <- read.delim("EPIC Cross sectional RP NEG Feature table.txt", skip=4, row.names = 1)
  
  #merge discriminant tables for heatmap rownames
  #for coffee
  #disctbl <- rbind(disctbl.pos, disctbl.neg)
  #saveRDS(disctbl, "Coffee discriminant table.rds")
  #disctbl <- readRDS("Coffee discriminant table.rds")
  
  #for alcohol
  #disctbl <- rbind(disc.alc.pos, disc.alc.neg)
  disctbl <- rbind(join.df.pos, join.df.neg)
  
  #need the same filter of peak table as for statistical analysis
  posfilt <- rowSums(ptpos > 1) > 265 #gives logical vector 
  negfilt <- rowSums(ptneg > 1) > 265
  posmat <- ptpos[ posfilt, ]
  negmat <- ptneg[ negfilt, ]
  
  #subset by discriminant feature number (transpose converts to matrix)
  discmatpos <- posmat[ disctbl.pos$feature,  ] %>% t
  discmatneg <- negmat[ disctbl.neg$feature,  ] %>% t

  discmat <- cbind(discmatpos, discmatneg)
  
  #filter the original peak table (for extraction of feature data)
  ptpos <- ptpos[disctbl.pos$feature, ]
  ptneg <- ptneg[disctbl.neg$feature, ]
  
  logmat   <- log2(discmat)
  colnames(logmat) <- paste(disctbl$mode, disctbl$med.int, disctbl$Mass.x, "@", round(disctbl$Rt, 3))
  #colnames(logmat) <- disctbl$feature
  cormat   <- cor(logmat)
  colnames(cormat) <- NULL
  
  require(corrplot)
  cplot <- corrplot(cormat, method="square", tl.col="black", order="hclust", tl.cex=0.45, 
           cl.ratio=0.2, cl.align="l", cl.pos="r", cl.cex=0.6, mar=c(1,1,1,1))
  
  disctbl$cluster <- as.numeric(rownames(cplot))
  
#---------------------------------------------------------------------------------------------------------
  
#manhattan plot updated for new pos and neg data.
#Load features with partial correlations
alc.all.ss <- readRDS("Alcohol features Manhattan.rds")
  
#Make colour vector for stripes, 6538 features
colvec <- c(rep(c(rep("col1", 500), rep("col2", 500)), 6), rep("col1", 538))
  
#Make data frame for Manhattan ggplot. Added variables: order, the factor of stripe colours, 
#conditional mutate of this for the final colour vector
df <- alc.all.ss %>% arrange(RT) %>% 
   mutate(order = 1:6538, 
        pointcol = as.factor(colvec), 
        direction = ifelse(Pcor > 0, "pos", "neg"), 
        absPcor = abs(Pcor),
        signif = ifelse(pval < 0.001 & direction == "pos", "col3", pointcol),
        signif2 = ifelse(pval < 0.001 & direction == "neg", "col4", signif))
           
#Call to ggplot (customise as needed)
ggplot(df, aes(x = order, y = absPcor, colour = signif2)) + geom_point() + theme_bw() +
  scale_colour_manual(values = c("darkgrey", "black", "red", "blue")) +
  scale_x_continuous(name = "Order of elution of spectral features (increasing lipophilicity)", 
                    expand = c(0.02,0)) +
  scale_y_continuous(name = "Partial Pearson correlation coefficient", expand = c(0.01,0.01)) +
  #geom_hline(yintercept = -0.00, colour = "white", size = 2) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
