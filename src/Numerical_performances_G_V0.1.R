#Calculate performance indexes: CC, AUC, PWSD,...

# Clean up everything
rm (list=ls())
setwd("/home/marco/workspace/CAGI_2015/p16_assessment/to_git/")

library(ROCR)
classThreshold <- c(65, 75, 90)

realData <- read.table(file= "data/step04/real_proliferation.txt", sep="\t", quote="", header=FALSE)
realData <- as.matrix(realData)
realValue <- as.numeric(realData[2:11, 3])


# Create empty results matrices
resultMatrix <- matrix(ncol = 5, nrow = 22)
colnames(resultMatrix) <- c("RMSE", "Pearson CC", "Pearson CC p-value", "Kendall CC", "Kendall CC p-value")
rownames(resultMatrix) <- rep("empty",22)
#
performanceMatrix <- matrix(ncol = 12, nrow = 22)
perfParameters <- c("Sens", "Spec","Bal. Accuracy","AUC" ) #,"#predictions_in_10_range"
listnames <-c()
for (l in classThreshold){
  for(m in perfParameters){
    listnames <-c(listnames, paste(m,l))
  }
}
colnames(performanceMatrix) <-listnames 
rownames(performanceMatrix) <-rep("empty",22)

#Calculate performance indexes: CC, AUC, PWSD,...
#load the prediction
for(i in 1:22) {
  prediction <- read.table(file=paste("data/step04/p16_Submission_", i , ".txt", sep=""), sep="\t", quote="", header=FALSE)
  prediction <- as.matrix(prediction)
  predictedValue <- as.numeric(prediction[2:11, 3]) * 100
  predictedSTD <- as.numeric(prediction[2:11, 4]) * 100 #as consequence, all "*" becomes "NA" NOTE: all Warnings messages are harmless

  
  #get the correlation coefficients
  if (T){
    pearsonCor.coeff <- cor(realValue , predictedValue , method = "pearson")
    pearsonCor.pval <- cor.test(realValue , predictedValue , method = "pearson")
    
    kendallCor.coeff <- cor(realValue , predictedValue , method = "kendall", use="pairwise")
    kendallCor.pval <- cor.test(realValue , predictedValue , method = "kendall", use="pairwise")
  }
  #get RMSE
  if (T){
    rmsd <- 0
    for(k in 1:length(realValue)) {
      rmsd <- rmsd + (predictedValue[k] - realValue[k])^2
    }
    rmsd <- sqrt(rmsd/length(realValue))
  }
  #print(paste("Submission", i, " rmsd: ", rmsd, " pearson: ",pearsonCor.coeff, " p-value: ", pearsonCor.pval[3], " kendall: ",kendallCor.coeff, " p-value: ", kendallCor.pval[3]))
  
  rownames(resultMatrix)[i] <- paste("Submission",i) 
  resultMatrix[i,] <-c(rmsd,pearsonCor.coeff, as.numeric(pearsonCor.pval[3]),kendallCor.coeff,as.numeric(kendallCor.pval[3]))
  
  #############################################
  #get AUC and other performance indexes
  if (T) {
    legendText <- as.character(classThreshold)
    
    #the set of AUC for every predictor
    performanceList <- c()
    #draw a ROC curve for every threshold
    for(j in 1:length(classThreshold)) {
      realClass <- realValue
      realClass[ which(realValue[] > classThreshold[j]) ] <- 1
      realClass[ which(realValue[] <= classThreshold[j]) ] <- 0
      
      #get predicted classes
      predClass <- predictedValue
      predClass[ which(predictedValue[] > classThreshold[j]) ] <- 1
      predClass[ which(predictedValue[] <= classThreshold[j]) ] <- 0
      #get pathogenic mutations
      positive <- which(realClass[] == 1)
      TP <- length( which(predClass[positive] == realClass[positive] ) )
      FN <- length( which(predClass[positive] != realClass[positive] ) )
      #get neutral mutations
      neutral <- which(realClass[] == 0)
      TN <- length( which(predClass[neutral] == realClass[neutral] ) )
      FP <- length( which(predClass[neutral] != realClass[neutral] ) )
      #get specificity and sensitivity
      sens <- TP / (TP + FN)
      spec <- TN / (TN + FP)
      bAccu <- (sens + spec) /2
      performanceList <- c(performanceList, sens, spec, bAccu)
      pred <- prediction( predictions= predictedValue, labels= realClass)
      perf <- performance(pred, "tpr", "fpr")
      auc <- performance(pred, "auc")
      auc <- unlist(slot(auc, "y.values"))
      legendText[j] <- paste(legendText[j], " (AUC: ", signif(auc,digits=3), ")")
      performanceList <- c(performanceList, auc)
      
    }
    ################################
    rownames(performanceMatrix)[i] <- paste("Submission",i) 
    performanceMatrix[i,] <- performanceList

    
  }
}

resultMatrix
performanceMatrix
finalMatrix <- cbind (resultMatrix,performanceMatrix)
finalMatrix <- round(finalMatrix,2)
finalMatrix

write.table(finalMatrix, "results/Tab3_SupTab3.tsv", sep = "\t" , quote = F, col.names = NA)
