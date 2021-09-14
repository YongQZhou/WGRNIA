#' Visualization of network
#' @description Show the gene regulator network which you provided. Visualization of gene regulatory network
#' @param grn A matrix which is the weighted adjacency matrix of the inferred network by this algorithm.
#' @param num The number of gene regulator imformation. Default: NULL.
#' @param file A connection, or a character string naming the file to write to. The default not save. Default: NULL.
#' @param TF The Transcription Factors (TFs) The default value NULL means that all the genes are used as candidate regulators. Default: NULL.
#' @param gene A vector of characters containing the target genes.
#' @details A labelled undirected network is plotted.
#' @return A D3 JavaScript force directed network graph.
#' @export
#'
#' @examples
#' B <- matrix(rnorm(1000,0,1),100)
#' TF <- paste0('G',1:20)
#' gene <- paste0('G',1:100)
#' grn <- run_aracne.a(B,gene,TF)
#' net <- network(grn,TF=TF,gene=gene)
network <- function(grn,num = NULL,TF,gene,file=NULL){
  grn <- as.matrix(grn)
  grn2 <- methods::as(grn,'sparseMatrix')
  row_pos <- grn2@i+1
  col_pos <- base::findInterval(seq(grn2@x)-1,grn2@p[-1])+1
  x <- grn2@x
  grn2 <- data.frame(i=row_pos,j=col_pos,x = x)
  if (is.null(num)) {
    num <- nrow(grn2)
  }
  grn2 <- grn2[order(abs(grn2$x),decreasing = T),]
  grn2 <- grn2[1:num,]
  gene2 <- gene[grn2$i]
  tf2 <- TF[grn2$j]
  #grn3 <- data.frame(TF = tf2, Gene = gene2, weight = grn2$x)

  name <- unique(c(TF,gene))
  TFid <- which(TF %in% name)
  geneid <- c(1:length(name))[-TFid]
  group <- c(1:length(name))
  group[TFid] <- 'TF'
  group[geneid] <- 'Gene'
  size <- c(1:length(name))
  size[TFid] <- 30
  size[geneid] <- 10
  grn4 <- data.frame(TF = grn2$j, Gene = grn2$i, weight = grn2$x)
  grn4[,1:2] <- grn4[,1:2]-1
  grn4$weight <- abs(grn4$weight) * 2
  datNodes <- data.frame(name = unique(c(TF,gene)),group = group,size = size)
  Colors <- 'd3.scaleOrdinal().range([ "#32CD32","#FFC0CB", "#8A2BE2", "#FF6347"])'
  net <- networkD3::forceNetwork(Links = grn4,Nodes = datNodes,Source = "TF",Target = "Gene",Value = "weight",
                                 NodeID = "name",Group = "group",Nodesize = "size" ,
                                 fontFamily="Times New Roman",fontSize = 20,linkColour="blue",
                                 charge = -100,opacity = 0.9,legend=T,arrows=T,
                                 colourScale = Colors,
                                 bounded=F,opacityNoHover=1.0,zoom = T)
  if (! is.null(file)) {
    networkD3::saveNetwork(net, file = file, selfcontained = TRUE)
  }
  return(net)
}

