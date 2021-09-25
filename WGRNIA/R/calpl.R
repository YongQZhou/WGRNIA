
#' Calculate the predictions of grn and the labels of the high-throughput golden standard.
#'
#' @description Evalute the infer network and high-throughput golden standard matrix. Calculate the predictions of grn and the labels of the high-throughput golden standard.
#' @param HGSmat A matrix which is the high-throughput golden standard. Rows contain genes and columns contain tfs ,the value between in one or zero.
#' @param grn A matrix which is the weighted adjacency matrix of the inferred network by this algorithm.
#' @param sym If set to TRUE, only the regulatory relationship in the high-throughput gold standard is considered, and the value of 0 is not considered, so that the number of true positive and false positive is equal, and both false negative and true negative are 0.
#' If set to FALSE, It is assumed that the regulatory relationships that do not exist in the high-throughput gold standard are not regulatory.The value as zero. Default: FALSE.
#' @param num The number of gene regulator imformation. Default: NULL
#'
#' @return A data.frame. The first element, dat$predictions, is a vector of numerical predictions. The second element, dat$labels, is a vector of cordatponding class labels.
#' @export
#'
#' @examples
#' A <- matrix(rnorm(200,0,1),20,10)
#' B <- matrix(rnorm(1000,0,1),100)
#' X <- matrix(0,100,20)
#' s <- c(1:10)
#' TF <- paste0('G',1:20)
#' gene <- paste0('G',1:100)
#' res <- demo_reg(A, B, s, X, 'ADMM',max.steps = 200, TF, gene, cl.core = 2,file=NULL,verbose=FALSE)
#' grn <- res$grn
#' HGSmat <- matrix(0,100,20)
#' HGSmat[sample(2000,200)] <- 1
#' dat <- calpl(HGSmat,grn)
calpl <- function(HGSmat,grn,sym = TRUE,num = NULL){
  grn <- abs(grn)
  HGSmat <- abs(HGSmat)
  grn <- as.matrix(grn)
  grn <- grn/max(grn)
  if (!is.null(num)) {
    spm1 = methods::as(grn, "CsparseMatrix")
    #  spm1 <- methods::as(grn,'sparseMatrix')
    row_pos <- spm1@i+1
    col_pos <- base::findInterval(seq(spm1@x)-1,spm1@p[-1])+1
    x <- spm1@x
    sm2 <- data.frame(i=row_pos,j=col_pos,x = x)
    
    sm2 <- sm2[order(sm2$x, decreasing = T), ]
    if (nrow(sm2) < num) { 
      sm2 <- sm2
    }else { sm2 <- sm2[1:num, ] }
    sm2 <- as.data.frame(sm2)
    grn2 <- Matrix::sparseMatrix(i = sm2$i, j = sm2$j, x = sm2$x,
                                 dims = c(nrow(grn), ncol(grn)))
    grn2 <- methods::as(grn2, "matrix")
    grn <- grn2
  }
  if (sym) {
    IND1 <- which(HGSmat == 1 & grn != 0, arr.ind = TRUE)
    N <- nrow(IND1)
    if (N != 0) {
      IND <- which(HGSmat == 0 & grn != 0, arr.ind = TRUE)
      N2 <- nrow(IND)
      if (N > N2) {
        N <- N2
      }
      set.seed(123)
      IND <- IND[sample(nrow(IND), N, replace = F), ]
      IND1 <- IND1[sample(nrow(IND1), N, replace = F),
      ]
      IND <- rbind(IND1, IND)
      dat <- grn[IND]
      dat <- matrix(c(dat, c(rep(1, N), rep(0, N))), ncol = 2)
      dat <- as.data.frame(dat)
      colnames(dat) <- c("predictions", "labels")
    }
    else {
      dat <- NULL
    }
  }
  if (!sym) {
    predictions <- as.numeric(grn)
    labels <- as.numeric(HGSmat)
    dat <- data.frame(predictions = as.numeric(predictions),
                      labels = as.numeric(labels))
  }
  return(dat)
}
