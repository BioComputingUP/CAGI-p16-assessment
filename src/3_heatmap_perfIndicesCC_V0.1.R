########################################################################################################################
#    this script calculate CC among performance indices and produce an heatmap figure: Figure 3 in assessment paper    #
########################################################################################################################
# Clean up everything
rm (list=ls())
setwd("/home/marco/workspace/CAGI_2015/p16_assessment/to_git/")

predictions <- read.delim("data/step05/predictions.txt", quote="")
predictions <- as.matrix(predictions)
colnames(predictions)[c(3,7,8)] <- c("RMSE","PWSD","PWSD10")

if (TRUE){
  #calculate the FULL correlation matrix
  matr <- matrix(0, ncol(predictions), ncol(predictions))
  for(i in 1:ncol(predictions)) {
    for(j in 1:ncol(predictions)) {
      matr[i, j] <- cor(predictions[,i] , predictions[,j] , method = "kendall", use="pairwise")
    }
  }
  #plot the correlation matrix
  require("plotrix")
  cellcolors<-matrix(NA, ncol(matr), ncol(matr))
  cellcolors[matr >= 0]<-
    color.scale(matr[matr >= 0],
                cs1=c(0.7,0),cs2=c(0,0.7),cs3=0)
  cellcolors[matr < 0]<-
    color.scale(matr[matr < 0],
                cs1=c(0.7,0),cs2=0,cs3=0.7)
  color2D.matplot(matr, show.values=TRUE, axes=FALSE, cellcolors=cellcolors, xlab="", ylab="", vcex = 1.6)
  axis(3, at=0.5:7.5,  labels=colnames(predictions),  cex.axis = 1.3, font = 2 ) 
  axis(2, at=0.5:7.5,  labels=rev(colnames(predictions)) , cex.axis = 1.3, font = 2 ) 
  #
  ppi = 600
  dev.copy(jpeg ,filename ="results/main_fig_3.jpeg",width=10*ppi, height=10*ppi, res=ppi)
  dev.off()
}
