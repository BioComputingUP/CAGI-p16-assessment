##################################################################################################################
#             This script calculate only PSWD10 and produce a table: Supp Table 2 in assessment paper            #
##################################################################################################################
# Clean up the environment
rm (list=ls())
setwd(".")

#get the real data
realData <- read.table(file="../data/real_proliferation.txt", sep="\t", quote="", header=FALSE)
realData <- as.matrix(realData)
realValue <- as.numeric(realData[2:11, 3])
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

diffTarget <- matrix(0,nrow = 23,ncol = 10)
#pswd10 <-c()
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
  #print(pval <0.05)
  #pswd10 <-c(pswd10,sum(pval > 0.05))
  diffTarget[i,] <- c(pval >0.05)
}

diffTarget[23,] <- apply(diffTarget,2, function(x) sum(x)) 
colnames(diffTarget) <- prediction[2:nrow(prediction),1]
diffTarget <- diffTarget[,c(4,5,6,7,8,1,2,3,9,10)]
rownames(diffTarget) <- c(paste("S", c(1:22)),"Total")
write.table(diffTarget, file='../results/sup_tab_2.txt', sep = "\t" , quote = F, col.names = NA )
