##################################################################################################################
#    this script calculate PCCs among predictions and produce an heatmap figure: Figure 2 in assessment paper    #
##################################################################################################################

setwd("/home/marco/workspace/CAGI_2015/p16_assessment/to_git/")
#get the real data
realData <- read.table(file = "data/step04/real_proliferation.txt", sep="\t", quote="", header=FALSE)
realData <- as.matrix(realData)
realValue <- as.numeric(realData[2:11, 3])
data <- realValue

#load the prediction
for(i in 1:22) {
  prediction <- read.table(file = paste("data/step04/p16_Submission_", i , ".txt", sep=""), sep="\t", quote="", header=FALSE)
  prediction <- as.matrix(prediction)
  predictedValue <- as.numeric(prediction[2:11, 3]) * 100
  data <- cbind(data, predictedValue)
}
#Fix colnames
colnames(data)[2:length(colnames(data))] <-   paste(colnames(data)[2:length(colnames(data))], c(1:22))
predictions <- data[,2:ncol(data)]
colnames(predictions) <- c(1:22) 

###############################################################
#draw heatmap of CC among all predictions  (ordered by group)
###############################################################

s1 <- c(10,16,21,22)
s2 <- c(2,11,17)
s3 <- c(5,12,18)
s4 <- c(8,14,19)
s5 <- c(9,15,20)
s6 <- c(6,13)
s7 <- c(3)
s8 <- c(7)
s9 <- c(4)
s10 <- c(1)
tot <-c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10) 
labID <- c("Yang&Zhou Lab.","Bromberg Lab.","BioFolD Lab.","Gough Lab.","Moult Lab.","Vihinen Lab.","Casadio Lab.","Dunbrack Lab.","Lichtarge Lab."," Anonymous.")
#####
matr_new <- matrix(0, ncol(predictions), ncol(predictions))
matr_new <- predictions[ , c(tot) ]

if (TRUE){
  matr <- matrix(0, ncol(matr_new), ncol(matr_new))
  for (i in 1:ncol(matr_new)) {  
    for (j in 1:ncol(matr_new)) {  
      matr[i, j] <- cor(matr_new[,i] , matr_new[,j] , method = "pearson", use="pairwise")
    }
  }
  matr <- round(matr, digits = 1)
  #plot the correlation matrix
  require("plotrix")
  par(mar=c(5.1, 4.1, 4.1,10), xpd=TRUE)
  par(pty="s")
  cellcolors<-matrix(NA, ncol(matr), ncol(matr))
  cellcolors[matr >= 0]<- color.scale(matr[matr >= 0],cs1=c(0.7,0),cs2=c(0,0.7),cs3=0)

  #CC between 0 and -1 will have the entire colour scale from red to dark blue
  if (TRUE){
    cellcolors[matr < 0]<- color.scale(matr[matr < 0], cs1=c(0,0.7),cs2=0,cs3=0.1)
  }

  color2D.matplot(matr, show.values=TRUE, axes=FALSE, cellcolors=cellcolors, xlab="", ylab="", vcex = 0.5)
  axis(3, at=0.5:21.5, labels = FALSE, cex.axis = 0.7, font = 1 ) 
  axis(2, at=0.5:21.5, labels = FALSE, cex.axis = 0.6, font = 1 ) 
  # set Tick Labels 
  if (TRUE){
    mtext(3, at=0.5:3.5,  text= s1,  cex.axis = 0.7, font = 1, col = "black", line =0.6)
    mtext(3, at=4.5:6.5,  text= s2,  cex.axis = 0.7, font = 1, col = "blue", line =0.6)
    mtext(3, at=7.5:9.5,  text= s3,  cex.axis = 0.7, font = 1, col = "magenta", line =0.6)
    mtext(3, at=10.5:12.5,  text= s4,  cex.axis = 0.7, font = 1, col = "red", line =0.6) 
    mtext(3, at=13.5:15.5,  text= s5,  cex.axis = 0.7, font = 1, col = "orange", line =0.6) 
    
    mtext(3, at=16.5:17.5,  text=s6,  cex.axis = 0.7, font = 1, col = "gold", line =0.6)
    mtext(3, at=18.5,  text= s7,  cex.axis = 0.7, font = 1, col = "pink", line =0.6) 
    mtext(3, at=19.5,  text= s8,  cex.axis = 0.7, font = 1, col = "yellow", line =0.6)
    mtext(3, at=20.5,  text= s9,  cex.axis = 0.7, font = 1, col = "green", line =0.6)
    mtext(3, at=21.5,  text= s10,  cex.axis = 0.7, font = 1, col = "grey", line =0.6)
  }
  if (TRUE){
    mtext(2, at=18.5:21.5,  text= rev(s1),  cex.axis = 0.7, font = 1, col = "black", line =0.6)
    mtext(2, at=15.5:17.5,  text= rev(s2),  cex.axis = 0.7, font = 1, col = "blue", line =0.6)
    mtext(2, at=12.5:14.5,  text= rev(s3),  cex.axis = 0.7, font = 1, col = "magenta", line =0.6)
    mtext(2, at=9.5:11.5,  text= rev(s4),  cex.axis = 0.7, font = 1, col = "red", line =0.6)
    mtext(2, at=6.5:8.5,  text= rev(s5),  cex.axis = 0.7, font = 1, col = "orange", line =0.6)

    mtext(2, at=4.5:5.5,  text=rev(s6),  cex.axis = 0.7, font = 1, col = "gold", line =0.6)
    mtext(2, at=3.5,  text= rev(s7),  cex.axis = 0.7, font = 1, col = "pink", line =0.6)
    mtext(2, at=2.5,  text= rev(s8),  cex.axis = 0.7, font = 1, col = "yellow", line =0.6)
    mtext(2, at=1.5,  text= rev(s9),  cex.axis = 0.7, font = 1, col = "green", line =0.6)
    mtext(2, at=0.5,  text= rev(s10),  cex.axis = 0.7, font = 1, col = "grey", line =0.6)
  }
  legend("topright", inset=c(-0.43,0.1), legend = labID, col=c("black", "blue", "magenta", "red", "orange", "gold", "pink", "yellow", "green", "grey"),
         pch = 18, cex=0.7, title = "Group identifiers")
  
  dev.copy(pdf, file = paste("results/main_fig_2.pdf", sep=""))
  dev.off()
}






