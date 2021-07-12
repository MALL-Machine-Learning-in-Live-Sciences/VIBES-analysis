draw_confusion_matrix <- function(cm) {
  require(viridis)
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  # create the matrix (Care with pos and  neg class, see in Confusion matrix the place ocupied)
  rect(150, 430, 240, 370, col=viridis(2)[1])
  text(195, 445, 'High', cex=1.2)
  rect(250, 430, 340, 370, col=viridis(3)[2])
  text(295, 445, 'Low', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=viridis(3)[2])
  rect(250, 305, 340, 365, col=viridis(2)[1])
  text(140, 400, 'High', cex=1.2, srt=90)
  text(140, 335, 'Low', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "Details", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1)
  text(30, 85, names(cm$byClass[2]), cex=1, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1)
  text(50, 85, names(cm$byClass[5]), cex=1, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1)
  text(70, 85, names(cm$byClass[6]), cex=1, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1)
  text(90, 85, names(cm$byClass[7]), cex=1, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1)
  text(70, 35, names(cm$overall[2]), cex=1, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1)
}  

