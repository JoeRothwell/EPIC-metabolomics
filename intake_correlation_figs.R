library(corrplot)
cplot <- corrplot(cormat, method = "square", tl.col = "black", order = "hclust", 
                  tl.cex = 0.5, cl.ratio = 0.2, cl.align = "l", cl.pos = "r", 
                  cl.cex = 0.6, mar = c(1,1,1,1))

# Get cluster order to paste into Excel
# writeClipboard(rownames(cplot))


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