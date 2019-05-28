library(tidyverse)
t1 <- readRDS("prepdata/OR_hcc_coffee_metabs.rds") %>% arrange(rowname)

#Generate sequence for row spacing
rowspace <- rev((1:32)[-seq(3, 32, by = 3)])

#Alternating 18 and 1 point styles
pointvec <- c(rep(c(18, 1), nrow(t1)/2))

dev.off()
#Sets margins (bottom, left, top, right)
par(mar=c(4,4,1,2))

library(metafor)
forest(slab = t1$rowname,
       ilab = t1$model, ilab.xpos = -1.5, ilab.pos = 4,
       x = t1$estimate, ci.ub = t1$conf.high, ci.lb = t1$conf.low, 
       refline = 1, rows = rowspace,
       xlab = "Odds ratio",
       ylim = c(1, max(rowspace) + 3), 
       xlim = c(-4, 6), 
       #alim = c(0.5, 2.5),
       pch = pointvec,
       #efac = 0,
       psize= rep(1.5, nrow(t1)),
       cex = 0.8)

text(6, max(rowspace) + 2, "OR [95% CI]", cex = 0.8, pos = 2)
text(c(-4, -1.5), max(rowspace) + 2, c("Compound", "Model"), cex = 0.8, pos=4)
