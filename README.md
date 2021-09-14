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
HGSmat[sample(ngene * ntf,0.2 * ngene * ntf)] <- 1
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
dat <- calpl(HGSmat = HGSmat,grn = grn,sym = T)
confu <- confusion(dat = dat)
print(confu)
## precision    recall  F1-score        TP        FP        FN 
##    0.5000    1.0000    0.6667   34.0000   34.0000    0.0000
auc <- calauc(dat = dat)
print(auc)
##  AUCROC   AUCPR 
## 0.59213 0.57890
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
##    algorithm precision recall F1-score  TP  FP FN  aucroc   aucpr       time
## 1    pearson       0.5      1   0.6667 104 104  0 0.47328 0.50122 0.52789402
## 2    pearson       0.5      1   0.6667 103 103  0 0.53747 0.57274 0.01795197
## 3    pearson       0.5      1   0.6667 156 156  0 0.49149 0.48740 0.02094388
## 4    pearson       0.5      1   0.6667 120 120  0 0.44837 0.45960 0.01695585
## 5    pearson       0.5      1   0.6667 144 144  0 0.52667 0.53229 0.01994491
## 6   aracne.a       0.5      1   0.6667  58  58  0 0.48915 0.48681 0.05485201
## 7   aracne.a       0.5      1   0.6667  55  55  0 0.61587 0.58341 0.02991891
## 8   aracne.a       0.5      1   0.6667  44  44  0 0.45868 0.47182 0.03094411
## 9   aracne.a       0.5      1   0.6667  63  63  0 0.48148 0.49732 0.02892303
## 10  aracne.a       0.5      1   0.6667  38  38  0 0.51004 0.55341 0.02988815
dat <- data.frame(alg = algorithm,auc = conlist$aucroc)
p <- draw_box(dat,'AUC')
print(p)
```

![img](https://github.com/YongQZhou/WGRNIA/blob/master/png/p3.png)





# Acknowledgements



# Contact

**Yongqiang Zhou**, postgraduate student at the School of Pharmaceutical Sciences (Shenzhen), Sun Yat-sen University.

Email: zhouyq67@mail2.sysu.edu.cn

