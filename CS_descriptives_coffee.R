# Coffee consumer stats. Prepare metadata from two datasets

library(tidyverse)
# Read in metadata from alcohol file (which has BMI data)
sumtab <- function() {
  alc_g <- read_csv("alcohol/alcohol.csv") %>% select(Idepic, Qe_Alc, R_BMI, Age_Dtq)
  
  #cs_metadata file now contains NAs where needed for easier subsetting. Read in and subset subjects (remove QC metadata)
  meta <- read.csv("data/cs_metadata.csv") %>% filter(present.pos == T & stype == "SA") %>%
    left_join(alc_g, by="Idepic") %>%  mutate(logQe_Alc = log(Qe_Alc + 1), 
      cupvols = case_when(country == "France" ~ 146.59, country == "Germany" ~ 209.32, 
                          country == "Greece" ~ 135.48, country == "Italy" ~ 55.2), cups = cof/cupvols)
  
  #Centre specific characteristics of the study population for paper
  meta <- meta %>% select(country, centre, sex, Age_Dtq, R_BMI, smoke, cof, cups)
  
  a <- meta %>% group_by(country, centre) %>% 
       summarise(n = n(), BMI = round(mean(R_BMI), 1), Cup_intake = round(median(cups),2),
       Cof_intake = round(median(cof)), Mean_age = round(mean(Age_Dtq), 1)) 
    
  b <- meta %>% group_by(centre, sex) %>% summarise(fem = n()) %>% filter(sex == "Female")
  c <- meta %>% group_by(centre, smoke) %>% summarise(smokers = n()) %>% filter(smoke == "Yes")
  d <- left_join(a, b, by = "centre")
  
  #Final summary table for paper
  sumtab <- left_join(d, c, by = "centre") %>%  select(Country=country, "EPIC Center"=centre, N=n, 
           "of which female" = fem, "Mean age" = Mean_age, 
           "Median coffee intake (mL)" = Cof_intake, "Mean BMI" = BMI, 
           "Current smokers at sample collection"=smokers)
}

tab <- sumtab()

#-------------------------------------------------------------------------------------------
# Visualisations
# histogram of all consumers
hist(meta$cof, breaks = 50, col = "gray")

# Boxplot of coffee intake by country (can change predictor to sex, smoke, alc. bev)
boxplot(sqrt(cof) ~ country, data=meta, col="dodgerblue", ylab="Alcohol intake (mL/day)")

# Boxplot for powerpoint slide. For publication reduce font size to 10 and B&W
library(ggplot2)
# ggplot(meta, aes(x=country, y=cups)) +
ggplot(meta, aes(x=country, y=cof)) + 
  #geom_boxplot(fill="grey", outlier.colour="white") +
  geom_boxplot(fill="dodgerblue", outlier.colour = "white") + 
  theme_bw(base_size = 10) + 
  geom_jitter(position=position_jitter(width=0.1), size=1) + #, colour="red")  +
  theme(panel.grid.minor=element_blank(), 
        #panel.grid.major = element_blank(),
        panel.border=element_rect(colour = "black")) +
  #ylim(c(0,1250)) +
  xlab("") + 
  ylab("Coffee intake mL/day") #+
  #ylab("Coffee intake cups/day")
  
ggsave("intake boxplot2.png", width=160, height = 120, units="mm")
ggsave("intake boxplot2.svg")

# For paper
ggsave("boxplot volume paper.png", width=80, height = 50, units="mm")
ggsave("boxplot cups paper.png", width=80, height = 50, units="mm")

# Boxplot for poster
library(latticeExtra)
bwplot(cof ~ country, data=meta, pch="|", ylab = "Coffee intake (mL/day)", ylim = c(-50, 1300),
       par.settings = list(box.rectangle = list(fill= "darkorange", col="black")))

# with tapply
cofsum1 <- tapply(meta$cof, meta$country, summary)
write.csv(cofsum, file="Intake summary by country.csv")

#-------------------------------------------------------------------------------------------------
# statistical tests
# Kruskal test (not as good as bootstrapped CIs as don't know which groups are different)
kruskal.test(meta$cof, meta$country)
kruskal.test(meta$cups, meta$country)
boxplot(meta$cups ~ meta$country)
fit <- lm(meta$cups ~ meta$country)
summary(fit)

# Medians and bootstrapped CIs of coffee intake or cups, 451 subjects
library(simpleboot)
tapply(meta$cof, meta$country, median)
tapply(meta$cups, meta$country, median)
b <- tapply(meta$cups, meta$country, FUN = one.boot, median, R=1000)
lapply(b, boot.ci, type=c("perc", "bca")) #two types output but percentile method used here

# Medians
# France  Germany   Greece    Italy 
# 248.3863 429.7500 140.0000  90.0000 

# Bootstrapped CIs
#          Percentile        BCa          
# France  (181.6, 320.0 )   (176.3, 310.8 ) 
# Germany (292.5, 450.0 )   (290.4, 450.0 )
# Greece  (112.0, 140.0 )   ( 70.5, 140.0 ) 
# Italy   (80.00, 92.31 )   (66.47, 90.00 ) 

# By base R with loop (warning, need to split data frame first, FR, D etc)
medians <- numeric(1000)
for (i in 1:1000) medians[i] <- median(sample(na.omit(FR$cof), replace = T))
ci <- quantile(medians, c(0.025, 0.975))
cat("95% confidence interval is (", ci, ")\n")

# ----

#Get coffee intake for all EPIC data
library(tidyverse)
library(haven)
foods <- read_csv("D:/dtqst_food_group.csv") %>% select(Country, Center, QGE130301, Idepic)

#Need to clean data after Python import of compressed SAS file
#See python script in MyScripts
#Practice code using data head only
hdfoods <- head(foods)
hdfoods %>% separate(Idepic, sep = "'", into = c("a", "b", "c"))

#On full data set
foods1 <- foods %>% separate(Idepic, sep = "'", into = c("a", "Idepic", "c")) %>% select(-a, -c)
class(foods1$Idepic)
saveRDS(foods1, "Coffee intakes all EPIC.rds")

intersect(foods1$Idepic, hcc$Idepic) %>% length()
dim(hcc)
