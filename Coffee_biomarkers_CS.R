# Coffee biomarker data from Profinder for CS study. Read in data generated from Profinder
# Samples, blanks and QCs were used (504 samples). "Repeat" files and 357 were omitted.

# Prepare data. Read in metadata and add cup volumes for countries

coffee_markers_cs <- function() {
  library(tidyverse)
  meta    <- read.csv("cs_metadata.csv")
  alc_g   <- read.csv("alcohol.csv") %>% select(Idepic, Qe_Alc, R_BMI)
  cupvols <- data.frame(country = levels(meta$country), cupvol = c(146.59, 209.32, 135.48, 55.2))
  meta    <- meta %>% left_join(cupvols, by="country") %>% mutate(cups = cof/cupvol)
  
  meta    <- meta %>% left_join(alc_g, by="Idepic")
  
  #Read in extra data for cyclo(prolyl-valyl)
  cyclo <- read_csv("CS cyclo pro val.csv") %>% slice(1)
  
  #read in pos and neg data for 10 biomarkers extracted in ProFinder. pos data:
  pos <- read_csv("CS 7 compounds pos.csv") %>% bind_rows(cyclo) %>% 
    select("Compound Name", contains("Area")) %>% t
  colnames(pos) <- pos[ 1, ]
  pos           <- pos[-1, ]
  
  pos.df <- cbind(ID = rownames(pos), tbl_df(pos)) %>% 
    separate(ID, into = c("Area", "datafile"), sep = " ", convert = T)
  
  #neg data
  neg <- read_csv("CS 3 compounds neg.csv") %>% select("Compound Name", contains("Area")) %>% t
  colnames(neg) <- neg[ 1, ]
  neg           <- neg[-1, ]
  
  neg.df <- cbind(ID = rownames(neg), tbl_df(neg)) %>% 
    separate(ID, into = c("Area", "datafile"), sep = " ", convert = T)
  
  #Join pos and neg data together, add metadata and filter subjects only
  
  posneg <- inner_join(pos.df, neg.df, by = "datafile") %>% 
    mutate(datalabel = paste(datafile, "(raw)", sep=""))
  
  wide <- inner_join(meta, posneg, by="datalabel") %>% filter(stype == "SA" & present.pos == T) %>%
    select(-Area.x, -Area.y)
  
  #Convert metabolite intensities to numeric. %<>% operator pipes and reassigns
  library(magrittr)
  wide[ , 37:47] %<>% lapply(function(x) as.numeric(x))
  
  #Make easier compound names (eg for STATA analysis)
  wide <- wide %>%
    rename( Trigonelline    = `Trigonelline (3-Carboxy-1-methylpyridinium betaine)`,
            #Quinic_acid     = `Quinic acid`, 
            #Hippuric_acid   = `Hippuric acid`, 
            #Cyclo_leu_pro   = `Cyclo(leucyl-prolyl)`,
            #Cyclo_isol_pro  = `Cyclo(isoleucyl-prolyl)`,
            #Cyclo_pro_val   = `Cyclo(prolyl-valyl)`,
            #Cat_sulf        = `Catechol sulfate`,
            AAMU            = `5-Acetylamino-6-amino-3-methyluracil (AAMU)`)
  
  #Remove hippuric acid, cyclo(leu-pro) and theophylline  
  #output <- wide %>% select(-`Cyclo_leu_pro`, -Hippuric_acid, -Theophylline)
}
wide <- coffee_markers_cs()

# Stats from extracted biomarkers

# get numbers of missing values
lapply(wide[, 37:47], function(x) sum(is.na(x))) #146 missing for quinic acid and 49 missing for AAMU

# logical indices for country subsetting
fr <- wide$country == "France"
d  <- wide$country == "Germany"
it <- wide$country == "Italy"
gr <- wide$country == "Greece"

#by(wide.df, wide.df$country, function(x) cor.test(wide.df$cof, wide.df[, 37:47]))

sum(fr)

# Calculate correlation coefficients with coffee
# for all countries
fit <- lapply(wide[  , 37:47], function(x) cor.test(wide[ , ]$cof, x, use = "pairwise.complete.obs"))
#lapply(wide.df[, 37:47], function(x) cor.test(x, wide.df$cups, use = "pairwise.complete.obs"))

# Caffeine/trigonelline ratio for AS
ratio <- wide %>% mutate(Caf_trig = Caffeine/Trigonelline)
cor.test(ratio$cof, ratio$Caf_trig, use = "pairwise.complete.obs")
plot(wide$Caffeine, wide$Trigonelline, col=wide$country, pch=18)

# or
library(psych)
fit <- corr.test(wide$cof, wide[, 37:47], use = "pairwise.complete.obs")


# Use broom to tabulate output of correlation
library(broom)
fit.summary <- lapply(fit, tidy)
#use do call to bind together the list of function arguments
cor.all <- do.call(rbind, fit.summary)

# multiple regression model using all biomarkers
cb <- wide %>% select(cof, Trigonelline:AAMU, -Cyclo_leu_pro, -Hippuric_acid, -Theophylline)

# get correlation matrix between biomarkers
bm.mat <- as.matrix(wide[, 37:47])
cormat <- cor(bm.mat, use = "pairwise.complete.obs")
cormat[cormat < 0] <- 0
colnames(cormat) <- rep(NA, 11)

library(corrplot)
library(RColorBrewer)
corrplot(cormat, method="color", col=colorRampPalette(c("blue","white","black"))(200),
         cl.lim=c(0,1), tl.col = "black", addgrid.col = "black")

# multiple regression model using all biomarkers
# maximum model
summary(lm(cb))
# minimum model (intercept only)
cb.lm <- lm(cof ~ 1, data = cb)
# add variables stepwise
add1(cb.lm, scope = cb)
cb.lm <- lm(cof ~ Trigonelline, data = wide)
cb.lm <- lm(cof ~ Trigonelline + Cyclo_pro_val + AAMU, data = cb)

fit  <- lm(cof ~ Trigonelline + Hippuric_acid + Theophylline + Paraxanthine + Caffeine +
             Cyclo_isol_pro + Cyclo_leu_pro + Cyclo_pro_val + Cat_sulf + AAMU, data = wide)


fit2 <- lm(cof ~ Trigonelline, data = wide)

summary(fit)
step(fit, direction = "both")
step(fit, direction = "backward")

fit.optimal <- lm(cof ~ Trigonelline + Hippuric_acid + Paraxanthine + Cyclo_leu_pro + Cyclo_pro_val + 
                    Cat_sulf + AAMU, data = wide)

library(MASS)
fit.rrr <- lm.ridge(cof ~ Trigonelline + Hippuric_acid + Theophylline + Paraxanthine + Caffeine +
           Cyclo_isol_pro + Cyclo_leu_pro + Cyclo_pro_val + Cat_sulf + AAMU, data = wide)

#-------------------------------------------------------------------------------------------------

# Analyse data by country. Convert to long, remove missings, add observation number, remove some outliers manually,
# remove Quinic acid (too many missings)

long.df <- wide %>% gather(Compound, intensity, -(order:datafile), na.rm = T) %>% 
  select(country, Compound, intensity, cof, cups) %>% mutate(int = intensity, obs = 1:n()) %>%
  filter(  obs != 2874 & obs != 1214 & obs != 906  & obs != 1665 & obs != 642  & obs != 554 & 
           obs != 1357 & obs != 728  & obs != 1158 & obs != 1012 & obs != 3018 & obs != 2939 & 
           obs != 2488 & obs != 2567 & obs != 2423 & obs != 4235 & obs != 3558,
           Compound != "Quinic acid")

# Faceted boxplots. Make facet labels
long.df$Compound <- factor(long.df$Compound, labels = c("AAMU", "Caffeine", "Catechol sulfate", 
  "Cyclo(isoleucyl-prolyl)", 
  "Cyclo(leucyl-prolyl)", 
  "Cyclo(prolyl-valyl)", 
  "Hippuric acid", 
  "Paraxanthine", 
  "Quinic acid", 
  "Theophylline", 
  "Trigonelline"))

ggplot(long.df, aes(x = country, y = int)) + 
  geom_boxplot(fill = "grey") +
  #geom_point(fill = "white") +
  #geom_errorbar(aes(ymax = Mn + StErr, ymin = Mn - StErr), width = 0.2) +
  xlab("Country") + ylab("Intensity (peak area)") +
  facet_wrap( ~ Compound, ncol = 2, scales = "free" ) +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), panel.border = element_rect(colour = "black"))

#---------------------------------------------------------------------------------------------
# scatter, all 451 points
ggplot(long.df, aes(x = cof, y = scale(int, center = F))) + 
  #ggplot(long.df, aes(x=cups, y=scale(intensity, center = F))) + 
  geom_point(aes(shape = country)) + theme_bw(base_size = 10) + 
  scale_shape_manual(values = c(0,1,3,4)) +
  facet_wrap( ~ Compound, ncol = 2, scales = "free_y") +
  geom_smooth(method = lm, se=F, colour="grey") +
  xlab("Median coffee intake, mL/day") + 
  #xlab("Median coffee intake, cups/day") +
  ylab("Scaled intensity") +
  theme(legend.position = "top", legend.title = element_blank())

ggsave("all points cups.png", width = 120, height = 180, units = "mm")
ggsave("all points vol.png", width = 120, height = 180, units = "mm")
#--------------------------------------------------------------------------------------------

# Median intensity vs median intake. Summarise data frame to get medians
sum.df <- long.df %>% group_by(country, Compound) %>% 
  summarise(MedInt = median(int), MedCof = median(cof), MedCups = median(cups),
            lowInt = quantile(int, 0.25),  highInt = quantile(int, 0.75),
            lowCof = quantile(cof, 0.25),  highCof = quantile(cof, 0.75),
            lowCup = quantile(cups, 0.25), highCup = quantile(cups, 0.75))

# Make facet labels for plot
sum.df$Compound <- factor(sum.df$Compound, labels = c("AAMU", "Caffeine", "Catechol sulfate", 
                                                        "Cyclo(isoleucyl-prolyl)", 
                                                        #"Cyclo(leucyl-prolyl)", 
                                                        "Cyclo(prolyl-valyl)", 
                                                        #"Hippuric acid", 
                                                        "Paraxanthine", 
                                                        "Quinic acid", 
                                                        #"Theophylline", 
                                                        "Trigonelline"))

# Plot (hash/unhash for volume or cups)
library(ggplot2)
#ggplot(sum.df, aes(x = MedCof, y = MedInt)) + 
ggplot(sum.df, aes(x = MedCups, y = MedInt)) + 
  geom_errorbar(aes(ymax = highInt, ymin = lowInt), colour="grey") +
  #geom_errorbarh(aes(xmax = highCof, xmin = lowCof), colour="grey") +
  geom_errorbarh(aes(xmax = highCup, xmin = lowCup), colour="grey") +
  geom_smooth(method = lm, se=FALSE, linetype="dotted", colour="black") +
  geom_point(aes(shape=country), size=2) + theme_bw() + 
  scale_shape_manual(values=c(15:18)) +
  facet_wrap( ~ Compound, ncol = 2, scales = "free_y") +
  #xlab("Median coffee intake / mL") + 
  xlab("Median coffee intake / cups") +
  ylab("Median intensity (counts)") +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.grid.major = element_blank() )

ggsave("med int vs med cups 8 cmpds.png", width=120, height = 180, units="mm")
ggsave("med int vs med vol 8 cmpds.png", width=120, height = 180, units="mm")

#---------------------------------------------------------------------------------------------



