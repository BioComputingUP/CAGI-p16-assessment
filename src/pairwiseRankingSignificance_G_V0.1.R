setwd("/home/marco/workspace/CAGI_2015/p16_assessment/to_git")
#print the pairwise significance of the P16 challenge evaluations

predictions <- read.delim("data/step05/predictionsRaw.txt", quote="")
predictions <- as.matrix(predictions)
                                                                     
for( i in c(9)) {  
  value <- predictions[, i]                                      
  finalRank <- rank( value, ties.method="first")
  
  #calculate the similarity of the ranking
  rankingSimilarity <- matrix(0, nrow(predictions), nrow(predictions))
  labelRankCorrelation <- vector(mode="character", length = nrow(predictions))
  for(k in 1:nrow(predictions)) {
    for(j in 1:nrow(predictions)) {
      tvar <- (value[k] - value[j]) / sqrt( 2*(sd(value)^2)/length(value) )
      rankingSimilarity[finalRank[k], finalRank[j]] <- 2 * pt(-abs(tvar),df=length(value) - 1)
      labelRankCorrelation[finalRank[k]] <- k
    }
  }

  
  #plot the ranking overlapping matrix
  require("plotrix")
  cellcolors<-matrix(NA, nrow(rankingSimilarity), nrow(rankingSimilarity))
  #get the upper triangular matrix
  index <- which(upper.tri(cellcolors, diag = FALSE))
  cellcolors[index] <- color.scale(rankingSimilarity[index], cs1=c(0, 1),cs2=c(0.7,1),cs3=c(0, 1))
  
  #the lower triangular matrix
  index <- which(lower.tri(cellcolors, diag = FALSE))
  cellcolors[index] <- color.scale(rankingSimilarity[index], cs1=c(0.7, 1),cs2=c(0,1),cs3=c(0, 1))
  
  color2D.matplot(rankingSimilarity, axes=FALSE, cellcolors=cellcolors, xlab="", ylab="", sub=paste("Significance of the ", colnames(predictions)[i]," ranking", sep=""))#show.values=TRUE, 
  axis(3, at=0.5:nrow(rankingSimilarity),  labels=labelRankCorrelation)
  axis(2, at=0.5:nrow(rankingSimilarity), labels=rev(labelRankCorrelation))
  dev.copy(jpeg, filename ="results/main_fig_4.jpeg ",  quality = 95, width = 800, height = 800)
  dev.off()
}


