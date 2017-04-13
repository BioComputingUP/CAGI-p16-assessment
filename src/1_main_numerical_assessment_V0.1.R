##############################################################################################################################
#                   This scripts calculate the main performance measures presented in main Tables 3 and 4                    #
#         it also writes 2 tabeles in /results that will be used by other scripts to produce the assessment figures          #
##############################################################################################################################

# Clean up the environment
rm (list=ls())
setwd(".")
library(ROCR)

#get the real data
realData <- read.table(file="../data/real_proliferation.txt", sep="\t", quote="", header=FALSE)
realData <- as.matrix(realData)
realValue <- as.numeric(realData[2:11, 3])


#create empty results matrices
resultMatrix <- matrix(ncol = 5, nrow = 22)
colnames(resultMatrix) <- c("RMSE", "Pearson CC", "Pearson CC p-value", "Kendall CC", "Kendall CC p-value")
rownames(resultMatrix) <- rep("empty",22)

performanceMatrix <- matrix(ncol = 12, nrow = 22)
perfParameters <- c("Sens", "Spec","Bal. Accuracy","AUC" )
# set 3 threshold for AUC
classThreshold <- c(65, 75, 90)

listnames <-c()
for (l in classThreshold){
  for(m in perfParameters){
    listnames <-c(listnames, paste(m,l))
  }
}
colnames(performanceMatrix) <-listnames 
rownames(performanceMatrix) <-rep("empty",22)

#############################################################
#     Calculate performance measures: CC, AUC, PWSD,...     #
#############################################################
#load the prediction
for(i in 1:22) {
  prediction <- read.table(file=paste("../data/p16_Submission_", i , ".txt", sep=""), sep="\t", quote="", header=FALSE)
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
    performanceList <- c()
    for(j in 1:length(classThreshold)) {
      realClass <- realValue
      realClass[ which(realValue[] > classThreshold[j]) ] <- 1
      realClass[ which(realValue[] <= classThreshold[j]) ] <- 0
      
      #get predicted classes
      classThreshold <- c(65, 75, 90)
      
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
      auc <- round(unlist(slot(auc, "y.values")), digits = 6)
      legendText[j] <- paste(legendText[j], " (AUC: ", signif(auc,digits=3), ")")
      performanceList <- c(performanceList, auc)
    }
    ################################
    rownames(performanceMatrix)[i] <- paste("Submission",i) 
    performanceMatrix[i,] <- performanceList
  }
}
#resultMatrix
#performanceMatrix
finalMatrix <- cbind (resultMatrix,performanceMatrix)
finalMatrix <- round(finalMatrix,2)
#finalMatrix
write.table(round(performanceMatrix[,c(1,2,3,5,6,7,9,10,11)], digits = 2), file='../results/sup_tab_3.txt', sep = "\t" , quote = F, col.names = NA )

############################################################
#                           PSWD                           #
############################################################
#the number of replicates for every prediction and experiment
nReplicates <- 10
realSTD <- as.numeric(realData[2:11, 4])

#the number of non * from the prediction file
count <- 0
#the cumulative std declared
cumStd <- vector(mode="numeric", length = 10)
for(i in 1:22) {
  prediction <- read.table(file=paste("../data/p16_Submission_", i , ".txt", sep=""), sep="\t", quote="", header=FALSE)
  prediction <- as.matrix(prediction)
  #check the first std in each prediction file
  if(prediction[2,4] != "*") {
    cumStd <- cumStd + as.numeric(prediction[2:11, 4]) * 100
    count <- count + 1
  }
}
cumStd <- cumStd / count

pswd <-c()
#load the prediction
for(i in 1:22) {
  prediction <- read.table(file=paste("../data/p16_Submission_", i , ".txt", sep=""), sep="\t", quote="", header=FALSE)
  prediction <- as.matrix(prediction)
  predictedValue <- as.numeric(prediction[2:11, 3]) * 100
  predictedSTD <- as.numeric(prediction[2:11, 4]) * 100
  predictedSTD[which(is.na(predictedSTD))] <- cumStd[which(is.na(predictedSTD))]
  predictedSTD[which(predictedSTD[] == 0)] <- cumStd[which(predictedSTD[] == 0)]
  #fix to 6 the number of replicates of each experiment
  tvar <- (realValue - predictedValue) / sqrt( (predictedSTD^2)/nReplicates + (realSTD^2)/nReplicates )
  pval <- pt(-abs(tvar),df=pmin(nReplicates, nReplicates)-1)
  pswd <-c(pswd,sum(pval > 0.05))
}

############################################################
#                         PSWD 10                          #
############################################################
pswd10 <-c()
#load the prediction
for(i in 1:22) {
  prediction <- read.table(file=paste("../data/p16_Submission_", i , ".txt", sep=""), sep="\t", quote="", header=FALSE)
  prediction <- as.matrix(prediction)
  predictedValue <- as.numeric(prediction[2:11, 3]) * 100
  predictedSTD <- as.numeric(prediction[2:11, 4]) * 100
  predictedSTD[] <- 10 
  #set to 1% STD when it is not defined or 0
  predictedSTD[which(is.na(predictedSTD))] <- cumStd[which(is.na(predictedSTD))]
  predictedSTD[which(predictedSTD[] == 0)] <- cumStd[which(predictedSTD[] == 0)]
  #fix to 6 the number of replicates of each experiment
  tvar <- (realValue - predictedValue) / sqrt( (predictedSTD^2)/nReplicates + (realSTD^2)/nReplicates )
  pval <- pt(-abs(tvar),df=pmin(nReplicates, nReplicates)-1)
  pswd10 <-c(pswd10,sum(pval > 0.05))
}
#############################################################

######              Compose assessment matrix          #####   
pRaw <- cbind(resultMatrix[,c(2,4,1)], performanceMatrix[,c(4,8,12)])
pRaw <- cbind(pRaw,pswd)
pRaw <- cbind(pRaw,pswd10)

######          Compose the final ranking matrix       #####
p <- apply(pRaw,2, function(x) rank(-x, ties.method = "min"))
p[,3] <- rank(pRaw[,3],ties.method = "min")
write.table(p, file='../results/predictions.txt', sep = "\t" , quote = F, row.names = F )

######       Compose the final assessment matrix      #####
pAverall <-p[,c(2,3,5,8)]
averall <- apply(pAverall,1,function(x) mean(x))
pRaw <-cbind(pRaw,averall)
write.table(pRaw, file='../results/predictionsRaw.txt', sep = "\t" , quote = F, row.names = F )

