#' draw_box
#' @description Show the plot of AUC or AUCPR by different algorithms of bootstrap.
#' @param dat The data.frame with the name of algorithm and the value of auc of aucpr.
#' @param name AUC or AUCPR
#'
#' @return The plot of AUC or AUCPR by different algorithm of bootstrap
#' @export
#' @import RColorBrewer ggplot2
#'
#' @examples
#' set.seed(123)
#' A <- matrix(rnorm(200,0,1),20,10)
#' B <- matrix(rnorm(1000,0,1),100)
#' X <- matrix(0,100,20)
#' s <- c(1:10)
#' TF <- paste0('G',1:20)
#' gene <- paste0('G',1:100)
#' HGSmat <- matrix(0,100,20)
#' HGSmat[sample(2000,500)] <- 1
#' algorithm <- c('pearson','spearman','bc3net','aracne.a','pcit','ADMM','OMP')
#' conlist <- grn_main(A=A,B=B,X=X,seq = s,gene=gene,TF=TF,algorithm = algorithm,
#'                     HGSmat = HGSmat,file=FALSE,verbose=TRUE)
#' alg = conlist$algorithm
#' auc = conlist$aucroc
#' dat <- data.frame(alg = alg,auc = auc)
#' p <- draw_box(dat,'AUC')
#' print(p)
draw_box <- function(dat,name = 'AUC'){
  colnames(dat) <- c('alg','auc')
  color2 <- RColorBrewer::brewer.pal(9,name = 'Set1')
  dat[,1] <- factor(dat[,1],levels = unique(dat[,1]))
  alg <- dat[,1]
  auc <- dat[,2]
  p <- ggplot(dat, aes(x = alg, y = auc, fill = alg))
  p2 <- p + theme_light() + geom_boxplot(width = 0.5,alpha=0.8,size=0.3,outlier.size = 0.5)  +
    labs(title = paste0("Plot of ",name," by different algorithm of bootstrap"), x = '', y = name,fill='methodtype')+
    scale_y_continuous(limits = c(0.3,1)) +
    theme(
      plot.title = element_text(hjust = 0.5,colour="black",face = 'bold',size = 15),
      axis.text.x=element_text(angle = 45,hjust = 0.5,vjust = 0.3, colour="black",size =13,face = 'bold'),
      axis.text.y=element_text(angle = 90,hjust = 1,vjust = 0.3, colour="black",size =13,face = 'bold'),
      axis.title.y=element_blank(),
      legend.title = element_text(colour="black",face = 'bold',size = 13),
      legend.text = element_text(colour="black",face = 'bold',size = 13),
      panel.border=element_rect(colour = "black")) +
    scale_fill_manual(values = unique(color2))
  p2 <- p2  + geom_hline(aes(yintercept = 0.5),colour='red',linetype = 'dashed')
  return(p2)
}
