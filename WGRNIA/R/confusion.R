
#' Confusion Matrix
#' @description Calculate the confusion martix, takes the weighted adjacency matrix of the inferred network by this algorithm and the true or high-throughput golden standard network for validation.
#' @param dat A data.frame. The first element, dat$predictions, is a vector of numerical predictions. The second element, dat$labels, is a vector of cordatponding class labels.
#'
#' @return A vector of confusion matrix including some ealuation parameters. (TP, FP, FN, precision, recall, F1-score).
#' @export
#'
#' @references G. Altay, F. Emmert-Streib, "Inferring the conservative causal core of gene regulatory networks", BMC Systems Biology, (2010) 4:132.
#' @examples
#' dat <- data.frame(predictions = runif(200, min=0, max=1),labels = c(rep(1,100),rep(0,100)))
#' score <- confusion(dat)
#'
confusion <- function(dat){
  predictions <- dat$predictions
  labels <- dat$labels
  predictions[predictions != 0] <- 1
  TP <- sum(predictions * labels)
  lab2 <- abs(labels - 1)
  FP <- sum(predictions * lab2)
  pred3 <- abs(predictions - 1)
  FN <- sum(pred3 * labels)
  #TN <- sum(pred3 * lab2)
  precision <- TP/(TP + FP)
  recall <- TP/(TP + FN)
  #specificity <- TN/(TN + FP)
  F1score <- 2 * TP/(2 * TP + FP + FN)
  output <- c(precision, recall, F1score, TP,
              FP, FN)
  namesv <- c("precision", "recall", "F1-score",
              "TP", "FP", "FN")
  names(output) <- namesv
  output <- round(output,4)
  return(output)

}
