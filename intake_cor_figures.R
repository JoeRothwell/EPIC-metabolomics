# Heatmap and Manhattan plot for correlations

# Correlation heatmap. First run intakecor_cs() function

discmat <- rbind(cofpos, cofneg) %>% select(starts_with("X")) %>% as.matrix %>% t
dat     <- rbind(cofpos, cofneg) %>% select(-starts_with("X"))
logmat  <- log2(discmat)
colnames(logmat) <- paste(dat$mode, dat$medint, dat$Mass, "@", round(dat$RT, 3))
#colnames(logmat) <- paste(disctbl$feature)
#colnames(logmat) <- disctbl$feature
cormat   <- cor(logmat)
colnames(cormat) <- NULL

library(corrplot)
cplot <- corrplot(cormat, method = "square", tl.col = "black", order = "hclust", 
                  tl.cex = 0.5, cl.ratio = 0.2, cl.align = "l", cl.pos = "r", 
                  cl.cex = 0.6, mar = c(1,1,1,1))

# Get cluster order to paste into Excel
# writeClipboard(rownames(cplot))

# Manhattan plot ----

manhattandata <- function() {
  library(tidyverse)
  library(RColorBrewer)
  cof <- readRDS("Coffee_features_Manhattan.rds")
  
  # Vector of colours for stripes, length 5934
  colvec <- rep(c("col1", "col2"), 6, each = 500)[1:nrow(cof)]
  colvec2 <- rep(brewer.pal(10, "Paired"), each = 600)[1:nrow(cof)]
  
  #Make data frame for Manhattan ggplot. Added variables: order, the factor of colours, 
  #conditional mutate of this for the final colour vector
  df <- cof %>% arrange(RT) %>% 
    mutate(order     = 1:n(), pointcol = colvec, 
           signif    = ifelse(pval < 0.05 & Pcor > 0, "col3", pointcol),
           signif2   = ifelse(pval < 0.05 & Pcor < 0, "col4", signif),
           signif3   = ifelse(pval < 0.05, "empty", "filled"),
           #mz2       = ifelse(pval > 0.9953192, paste("m/z", Mass, sep = " "), NA)
           mz        = ifelse(abs(Pcor) > 0.30, paste("m/z", Mass, sep = " "), NA))
}
df <- manhattandata()

#Calculate raw and adjusted cut points for plot
#signif <- filter(df, pval < 0.05)

#For publication, reduce point size and font size
library(ggplot2)
ggplot(df, aes(x = order, y = abs(Pcor), colour = signif2, shape = signif3)) + 
  #ggplot(df, aes(x = order, y = -log10(rawp), colour = signif2, shape = signif3)) +
  geom_hline(yintercept = c(0.2, 0.4), linetype = c("dashed"), colour = "grey") +
  #geom_hline(yintercept = c(-log10(0.05), -log10(0.9953192)), linetype = "dashed", colour = "black") +
  geom_point(size = 1.2) + theme_bw(base_size = 10) +
  scale_colour_manual(values = c("darkgrey", "black", "red", "blue")) +
  #scale_colour_manual(values = c("darkgrey", "black", "black", "black")) +
  scale_shape_manual(values = c(1, 16)) + 
  scale_x_continuous(name = "Elution order (increasing lipophilicity)", expand = c(0.02,0)) +
  scale_y_continuous(name = "Partial Pearson correlation coefficient", expand = c(0.01, 0.01)) +
  #scale_y_continuous(name = "-log10(pvalue)", expand = c(0.01, 0.01)) +
  geom_text(aes(label = mz), hjust = -0.1, vjust = 0, size = 1.5, colour = "black") +
  #geom_hline(yintercept = -0.00, colour = "white", size = 2) +
  theme(legend.position = "none", 
        text = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Plot for publication one-column
ggsave("manhattan coffee paper1.png", width = 100, height = 60, units = "mm")
ggsave("manhattan coffee paper1.svg", width = 100, height = 60, units = "mm")


### Manhattan plot for coffee updated for new pos and neg data ---------------------------------------
# From old script
#Load features with partial correlations
library(tidyverse)
disc.cof <- readRDS("Coffee_features_Manhattan.rds")

#Make colour vector for stripes, 5934 features
colvec   <- c(rep(c(rep("col1", 500), rep("col2", 500)), 5), rep("col1", 500), rep("col2", 434))

#Make data frame for Manhattan ggplot. Added variables: order, the factor of stripe colours, 
#conditional mutate of this for the final colour vector
df <- disc.cof %>% arrange(RT) %>% 
  mutate(order     = 1:5934, 
         pointcol  = as.factor(colvec), 
         direction = ifelse(Pcor > 0, "pos", "neg"), 
         absPcor   = abs(Pcor),
         signif    = ifelse(pval < 0.05 & direction == "pos", "col3", pointcol),
         signif2   = ifelse(pval < 0.05 & direction == "neg", "col4", signif),
         signif3   = ifelse(pval < 0.05, "empty", "filled"),
         mz        = ifelse(absPcor > 0.30, paste("m/z", Mass, sep = " "), NA),
         mz2       = ifelse(pval > 0.9953192, paste("m/z", Mass, sep = " "), NA))

#Calculate raw and adjusted cut points for plot
#signif <- filter(df, pval < 0.05)

#Call to ggplot (customise as needed)
#For publication, reduce point size and font size
library(ggplot2)
ggplot(df, aes(x = order, y = absPcor, colour = signif2)) + 
  #ggplot(df, aes(x = order, y = -log10(rawp), colour = signif2, shape = signif3)) +
  geom_hline(yintercept = c(0.2, 0.4), linetype = c("solid"), colour = "grey") +
  #geom_hline(yintercept = c(-log10(0.05), -log10(0.9953192)), linetype = "dashed", colour = "black") +
  geom_point(size = 1) + theme_bw(base_size = 10) +
  scale_colour_manual(values = c("darkgrey", "black", "red", "blue")) +
  #scale_colour_manual(values = c("darkgrey", "black", "black", "black")) +
  scale_shape_manual(values = c(1, 16)) + 
  scale_x_continuous(name = "Elution order (increasing lipophilicity)", expand = c(0.02,0)) +
  scale_y_continuous(name = "Partial Pearson correlation coefficient", expand = c(0.01, 0.01)) +
  #scale_y_continuous(name = "-log10(pvalue)", expand = c(0.01, 0.01)) +
  geom_text(aes(label = mz2), hjust = -0.1, vjust = 0, size = 1.5, colour = "black") +
  geom_hline(yintercept = -0.00, colour = "white", size = 2) +
  theme(legend.position = "none", 
        text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#Plot for publication one-column
ggsave("manhattan coffee paper1.png", width = 100, height = 60, units = "mm")
ggsave("manhattan coffee paper1.svg", width = 100, height = 60, units = "mm")


