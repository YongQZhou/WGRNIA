# WGRNIA

Wrapping gene regulatory network inference algorithms.

# Introduction

In this work, We have developed a comprehensive open source R package WGRNIA(Wrapping gene regulatory network inference algorithms), based on the real expression data and clearly defined benchmark dataset, aimed at helping users of gene regulatory network construction tools to evaluate the performance under experimental data. This package contains of several inference algorithms that cover models of Mutual Information, Correlation, Machine learning, Regression, Greedy algorithms, including microarray and single-cell RNA-seq data of mouse embryonic stem cells. This package provides modules that can be input into gold standards as well as algorithms, which is helpful for the rigorous and scalable evaluation of the GRN construction methods under the specific data.

# Installation

### Suggest packages

Some gene regulation network inference algorithms needs to be installed manually.

1. We import 'ennet' from github to run 'ennet', the gene regulation network  inference algorithms.

   The detail to install the 'ennet', see from https://github.com/slawekj/ennet

2. We import 'foba' from CRAN to run 'FoBa', the gene regulation network  inference algorithms.

   `devtools::install_url('https://cran.r-project.org/src/contrib/Archive/foba/foba_0.1.tar.gz')`

### Install WGRNIA

`install.packages('devtools')`

`devtools::install_github('YongQZhou/WGRNIA/WGRNIA')`



# Usage for WGRNIA

## Preparing data

Expression matrix B and gene list need to be prepared. A random set of data is generated to demonstrate how to use this package.

```
library(WGRNIA)
ntf <- 30
ngene <- 100
nsample <- 10
set.seed(123)
A <- matrix(rnorm(ntf * nsample,0,1),ntf,nsample)
TF <- paste0('G',1:ntf)
B <- matrix(rnorm(ngene * nsample,0,1),ngene,nsample)
X <- matrix(0,ngene,ntf)
s <- c(1:10)
gene <- paste0('G',1:ngene)
HGSmat <- matrix(0,ngene,ntf)
HGSmat[sample(ngene * ntf,0.1 * ngene * ntf)] <- 1
```

## Network inference

Select the network inference method.If the regression method is selected, additional matrix A and matrix X need to be provided, and other methods do not need this datas.

```
res <- demo_reg(A, B, s, X, 'ADMM',max.steps = 200,TF, gene, cl.cores = 2,file=NULL,verbose=TRUE)
## The 1 -th sparsity is finished.
## The 2 -th sparsity is finished.
## The 3 -th sparsity is finished.
## The 4 -th sparsity is finished.
## The 5 -th sparsity is finished.
## The 6 -th sparsity is finished.
## The 7 -th sparsity is finished.
## The 8 -th sparsity is finished.
## The 9 -th sparsity is finished.
## The 10 -th sparsity is finished.
## [1] "ADMM has completed"
grn <- res$grn 
grn[1:5,1:5]
##           [,1]      [,2]      [,3]  [,4]      [,5]
## [1,] 0.1111111 0.0000000 0.0000000 0.000 0.1428571
## [2,] 0.1250000 0.0000000 0.0000000 0.000 0.0000000
## [3,] 0.2000000 0.0000000 0.0000000 0.125 0.1666667
## [4,] 0.1111111 0.0000000 0.1428571 0.000 0.2500000
## [5,] 0.1000000 0.1428571 0.0000000 0.000 0.2000000
res <- demo_other(B,'pearson',TF=TF,gene=gene,file=NULL,verbose=TRUE)
## [1] "pearson has completed"
grn <- res$grn 
grn[1:5,1:5]
##           G1        G2 G3 G4 G5
## G1 0.0000000 0.7148843  0  0  0
## G2 0.7148843 0.0000000  0  0  0
## G3 0.0000000 0.0000000  0  0  0
## G4 0.0000000 0.0000000  0  0  0
## G5 0.0000000 0.0000000  0  0  0
```

## Evalution

Calculate evaluation parameters FP, FN, precision, recall, F1-score, AUC and AUCPR.

```
dat <- calpl(HGSmat = HGSmat,grn = grn,sym = F)
confu <- confusion(dat = dat)
print(confu)
## precision    recall  F1-score        TP        FP        FN 
##    0.1538    0.0733    0.0993   22.0000  121.0000  278.0000
auc <- calauc(dat = dat)
print(auc)
##  AUCROC   AUCPR 
## 0.51414 0.11190
```

## ROC and PR curve

```
p1 <- draw.plot(dat,method = 'roc',algname = 'pearson')
p2 <- draw.plot(dat,method = 'pr',algname = 'pearson')
print(p1)
```
![img](https://github.com/YongQZhou/WGRNIA/blob/master/png/p1.png)

```
print(p2)
```

![img](https://github.com/YongQZhou/WGRNIA/blob/master/png/p2.png)

## Network diagram

Generate a visual network diagram and generate a data format that can be used for Cytoscape to establish the network.

```
net <- network(grn,TF = TF,gene = gene)
#print(net)
grn <- out_grn(grn,TF = TF,gene = gene)
head(grn)
##    TF Gene     weight
## 1  G9  G92 -0.8970992
## 2  G9  G13 -0.8969077
## 3 G13   G9 -0.8969077
## 4  G7  G97  0.8825489
## 5 G10  G59 -0.8575120
## 6  G1  G92  0.8438612
```

# The pipeline for WGRNIA to evalute this methods

```
algorithm <- c('pearson','aracne.a','mrnet','pcit','ADMM')
conlist <- grn_main(A=A,B=B,X=X,seq = s,gene=gene,TF=TF,algorithm = algorithm,
                    HGSmat = HGSmat,file=TRUE, verbose=FALSE,bootstrap_num = 5)
print(conlist[1:10,])
##    algorithm precision recall F1-score TP  FP  FN  aucroc   aucpr       time
## 1    pearson    0.1207 0.1767   0.1434 53 386 247 0.51618 0.11255 0.12199903
## 2    pearson    0.1283 0.2133   0.1602 64 435 236 0.52880 0.12739 0.04601097
## 3    pearson    0.0952 0.2500   0.1379 75 713 225 0.49140 0.09997 0.06400490
## 4    pearson    0.0901 0.1767   0.1194 53 535 247 0.48994 0.09882 0.04604197
## 5    pearson    0.1149 0.2500   0.1574 75 578 225 0.51953 0.11320 0.06100011
## 6   aracne.a    0.1347 0.1100   0.1211 33 212 267 0.51580 0.11260 0.04399896
## 7   aracne.a    0.1565 0.1200   0.1358 36 194 264 0.52538 0.12646 0.03797507
## 8   aracne.a    0.1041 0.0767   0.0883 23 198 277 0.50128 0.10151 0.04399300
## 9   aracne.a    0.1048 0.1100   0.1073 33 282 267 0.50191 0.10274 0.04399896
## 10  aracne.a    0.0845 0.0600   0.0702 18 195 282 0.49385 0.09908 0.04300499
dat <- data.frame(alg = algorithm,auc = conlist$aucroc)
p <- draw_box(dat,'AUC')
print(p)
```

![img](https://github.com/YongQZhou/WGRNIA/blob/master/png/p3.png)



# Acknowledgements



# Contact

**Yongqiang Zhou**, postgraduate student at the School of Pharmaceutical Sciences (Shenzhen), Sun Yat-sen University.

Email: zhouyq67@mail2.sysu.edu.cn

