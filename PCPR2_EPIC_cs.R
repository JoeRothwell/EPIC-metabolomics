# Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
library(MetabolAnalyze)
library(car)

set.seed(111)
mat <- matrix(rnorm(10000), nrow = 100)
meta <- data_frame(var1 = rnorm(1:100), var2 = rnorm(1:100), var3 = rnorm(1:100))

source("Data prep.R")
#selected <- meta[meta$present.pos == T, ]
mat <- pt %>% as.matrix
mat <- ifelse(mat == 0 | is.na(mat), 1, mat)
X_MatScaled <- mat #X matrix
Z_Meta <- selected[, 17:19] #Y variables
Z_Meta <- meta

# Center the data / Scale the data, edit the parameter "pareto" or "unit"

X_MatCentered <- scale(X_Mat, center = TRUE, scale = FALSE) 
X_MatScaled <- scaling(X_MatCentered, type = "pareto")

Z_MetaRowN <- nrow(Z_Meta)
Z_MetaColN <- ncol(Z_Meta) 
ColNames <- names(Z_Meta)

# Perform PCA and retain no. components that account for threshold
# Obtain eigenvectors

pca <- prcomp(X_MatScaled) 
  
pct_threshold <- .8 # Amount of variability desired to be explained

X_MatScaled_t <- t(X_MatScaled)

Mat2        <- X_MatScaled %*% X_MatScaled_t
eigenData   <- eigen(Mat2)
eigenValues <- eigenData$values
ev_n        <- length(eigenValues)

eigenVectorsMatrix <- eigenData$vectors
percents_PCs <- eigenValues / sum(eigenValues)
my_counter_2 <- 0
my_sum_2 <- 1

for (i in ev_n:1) {
  my_sum_2 <- my_sum_2 - percents_PCs[i]
  if (my_sum_2 <= pct_threshold ) my_counter_2 <- my_counter_2 + 1
}

if(my_counter_2 < 3) pc_n = 3 else pc_n = my_counter_2

pc_data_matrix <- matrix(data = 0, nrow = Z_MetaRowN * pc_n, ncol = 1) 
mycounter <- 0

for (i in 1 : pc_n){
  for (j in 1 : Z_MetaRowN){ 
    mycounter <- mycounter + 1                          
    pc_data_matrix[mycounter, 1] <- eigenVectorsMatrix[j, i] 
  }
}

AAA <- Z_Meta[rep(1 : Z_MetaRowN, pc_n), ]
Data <- cbind(pc_data_matrix, AAA)

#Perform linear multiple regression models on each eigenvector with factors of interest as explanatory variables
#Categorical variables should be processed by as.factor
#To be edited with your factors names

#Z_Meta$Altitude <- Z_Meta$Altitude
DataCol <- ncol(Data)

type3mat <- matrix(data = 0, nrow = pc_n, ncol = DataCol ) 

ST_ResidualR2 <- matrix(data = 0, nrow = pc_n, ncol = 2)

for (i in 1:pc_n){ 
  yy <- (i-1) * Z_MetaRowN 
  y  <- yy + 1                                                                                                                      
  TotSumSq <- var(Data[y : yy + Z_MetaRowN, 1]) * (Z_MetaRowN - 1)
  
  #Edit the linear model with your factors
  Model            <- lm(pc_data_matrix ~  var1 + var2 + var3, 
                         Data[y : yy + Z_MetaRowN, ])
  AnalysisVariance <- Anova(Model, type = 3)
  SumSq            <- AnalysisVariance[1] 
  Residuals        <- SumSq[DataCol + 1, ] 
  RR               <- Residuals/TotSumSq
  R2               <- 1 - RR
  ST_ResidualR2[i, ]  <- c(R2, RR)
  colnames(ST_ResidualR2) <- c("ST_R2", "ST_Residuals")
  
  for (j in 1 : DataCol){
    type3mat[i,j] <- as.numeric(SumSq[j + 1, 1])
    colnames(type3mat) <- c(ColNames, "SumSqResiduals")
  } 
  
}

partialR2Matrix <- matrix(data = 0, nrow = pc_n, ncol = DataCol-1 )

for (i in 1:pc_n){
  for (j in 1:(DataCol - 1)) partialR2Matrix[i,j] = type3mat[i,j] / (type3mat[i, DataCol] + type3mat[i,j]) 
}

partialR2MatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = DataCol)  

for (i in 1:pc_n){
  weight <- eigenValues[i]/sum(eigenValues[1 : pc_n]) 
  for (j in 1:DataCol - 1){
    partialR2MatrixWtProp[i, j      ] <- partialR2Matrix[i, j] * weight
    partialR2MatrixWtProp[i, DataCol] <- ST_ResidualR2[i,1]*weight 
  }
}

pR2Sums <- colSums(partialR2MatrixWtProp) * 100 
plotnames <- c( ColNames, "R2")
bp        <- barplot(pR2Sums, xlab = "Factors", ylab = "Weighted Rpartial2", ylim= c(0,60),col = c("red"), las=2)
axis(1, at = bp, labels = plotnames, xlab = "Factors", cex.axis = 0.5, las=2) 
values     <- pR2Sums
new_values <- round(values, 3)
text(bp, pR2Sums, labels = new_values, pos=3, cex = 0.8)

output <- data.frame(plotnames, pR2Sums)

#--------------------------------------------------------------------------------------------------

library(ggplot2)
ggplot(data=output, aes(x=plotnames, y=pR2Sums)) + geom_bar(stat="identity", colour="black", fill="red") +
  theme_bw() + geom_text(aes(label=new_values),  vjust=-0.5) + ylim(0,60) +
  scale_x_discrete(name="Factor", limits=c("Replicate", "Experiment", "Altitude", "R2"))

#library(lattice)
#barchart(pR2Sums~plotnames, data=output, origin=0, col="pink")