#' Demo other - The demo of other method
#' @description This is the main function to call different algorithms for network inference.
#'
#' @param datexpr Expression matrix (genes x samples). Every row is a gene, every column is a sample.
#' @param algorithm The algorithm of inferring grn. We provide 13 algorithms in this function including 'pearson','spearman','kendall','WGCNA','bc3net','c3net','aracne.a','aracne.m',
#' 'clr','mrnet','mrnetb','GeneNet','Genie3.ET','Genie3.RF','pcit'. Default: pearson.
#' @param TF The Transcription Factors (TFs) The default value NULL means that all the genes are used as candidate regulators. Default: NULL.
#' @param gene A vector of characters containing the target genes.
#' @param cl.cores The number of cores in computer that you want to utilize. Default: 2.
#' @param file A connection, or a character string naming the file to write to. The default not save the grn matrix. Default: NULL.
#' @param verbose If set to TRUE, a feedback on the progress of the progress of the calculations is given. Default: TRUE.
#' @return A list with the function of parameters ,algorithms, elapse time, grn which a matrix which is the weighted adjacency matrix of the network inferred by the this algorithm
#' @export
#' @details The demo function aims to call algorithms to infer gene regulatory network.
#' @examples
#' datexpr <- matrix(rnorm(1000,0,1),100)
#' TF <- paste0('G',1:20)
#' gene <- paste0('G',1:100)
#' res <- demo_other(datexpr=datexpr,algorithm = 'pearson' ,TF=TF,gene=gene,file=NULL,verbose=TRUE)
#' grn <- res$grn
demo_other <- function(datexpr, algorithm='pearson', TF=NULL, gene, cl.cores = 2,file = NULL,verbose = TRUE){
  call <- match.call()
  path1 <- getwd()
  if (! is.null(file)) {
    dir <- 'Output'
    if (! dir.exists(dir)) {
      dir.create(dir)
    }
    Outpath <- str_c(dir,'/Output_', algorithm)
    if (! dir.exists(Outpath)) {
      dir.create(Outpath)
    }
    setwd(Outpath)
  }

  colnames(datexpr) <- NULL
  rownames(datexpr) <- gene

  if (algorithm == 'pearson') {
    t1 <- Sys.time()
    grn =  cor_pearson(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }

  if (algorithm == 'spearman') {
    t1 <- Sys.time()
    grn = cor_spearman(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }

  if (algorithm == 'kendall') {
    t1 <- Sys.time()
    grn = cor_kendall(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }


  if (algorithm == 'WGCNA') {
    t1 <- Sys.time()
    grn = run_wgcna(datexpr,gene,TF,file,cl.cores = cl.cores)
    t2 <- Sys.time()
  }

  if (algorithm == 'bc3net') {
    t1 <- Sys.time()
    grn = run_bc3net(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }

  if (algorithm == 'c3net') {
    t1 <- Sys.time()
    grn = run_c3net(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }

  if (algorithm == 'aracne.a') {
    t1 <- Sys.time()
    grn = run_aracne.a(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }

  if (algorithm == 'aracne.m') {
    t1 <- Sys.time()
    grn = run_aracne.m(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }

  if (algorithm == 'clr') {
    t1 <- Sys.time()
    grn = run_clr(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }


  if (algorithm == 'mrnet') {
    t1 <- Sys.time()
    grn = run_mrnet(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }


  if (algorithm == 'mrnetb') {
    t1 <- Sys.time()
    grn = run_mrnet(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }

  if (algorithm == 'GeneNet') {
    t1 <- Sys.time()
    grn = run_GeneNet(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }

  if (algorithm == 'Genie3.ET') {
    t1 <- Sys.time()
    grn = run_GENIE3.ET(datexpr,gene,TF,file,cl.cores)
    t2 <- Sys.time()
  }

  if (algorithm == 'Genie3.RF') {
    t1 <- Sys.time()
    grn = run_GENIE3.RF(datexpr,gene,TF,file,cl.cores)
    t2 <- Sys.time()
  }

  if (algorithm == 'pcit') {
    t1 <- Sys.time()
    grn = run_pcit(datexpr,gene,TF,file)
    t2 <- Sys.time()
  }


  if (algorithm == 'ennet') {
    t1 <- Sys.time()
    grn = run_ennet(datexpr,gene,TF,file,cl.cores = cl.cores)
    t2 <- Sys.time()
  }

  if (verbose) {
    print(sprintf("%s has completed", algorithm))
  }
  time <- as.numeric(t2) - as.numeric(t1)
  if (is.null(file)) {
    object <- list(call = call, algorithm = algorithm, time = time, grn = grn)
  }else{
    object <- list(call = call, algorithm = algorithm, time = time, savePath = Outpath,grn = grn)
  }

  setwd(path1)
  return(object)

}


