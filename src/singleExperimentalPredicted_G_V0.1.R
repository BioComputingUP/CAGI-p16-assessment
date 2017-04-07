#the script draw the experimental vs predicted graph, which compares the value of each prediction and the value of the corresponding real data

setwd("/home/marco/workspace/CAGI_2015/p16_assessment/to_git/")
workingFolder <- "data/step04/"
#get the real data
realData <- read.table(file = "data/step04/real_proliferation.txt", sep="\t", quote="", header=FALSE)
realData <- as.matrix(realData)
realValue <- as.numeric(realData[2:11, 3])

if (TRUE){
  #load the prediction
  for(i in 1:22) {
    prediction <- read.table(file=paste("data/step04/p16_Submission_", i , ".txt", sep=""), sep="\t", quote="", header=FALSE)
    prediction <- as.matrix(prediction)
    predictedValue <- as.numeric(prediction[2:11, 3]) * 100
    #reset boundaries
    predictedValue[which(predictedValue < 50)] <- 50
    predictedValue[which(predictedValue > 100)] <- 100
    par(srt = 45, mar=c(5,5,2,2))
    plot(realValue, predictedValue, main=  paste("Submission", i), pch = 20, col = "#0000FFFF", xlim=c(50, 103), ylim=c(50, 103), xlab="Experimental values", ylab="Predicted values", cex.axis=1.3, cex.lab=1.5, cex.main=1.5, cex = 2)
    text(realValue, predictedValue, labels = realData[2:11, 1], col = "#0000FFFF", pos = 4, cex=1.3)
    segments(48, 48, 105, 105, col= "#0000FFFF", lty= 2, lwd = 2)
    dev.copy(pdf, file = paste( "results/ExperimentalVSPredicted_", i ,".pdf", sep="")) #, width = 60, height = 60
    dev.off()
  }
}
