## [Machine Learning–Based Differential Network Analysis: A Study of Stress-Responsive Transcriptomes in **Arabidopsis**](https://github.com/cma2015/mlDNA)

## Introduction

mlDNA is an machine learning (ML) based methodology for transcriptome analysis via comparison of gene coexpression networks. The mlDNA first used a ML-based filtering process to remove nonexpressed, constitutively expressed, or non-stress-responsive “noninformative” genes prior to network construction, through learning the patterns of 32 expression characteristics of known stress-related genes. The retained “informative” genes were subsequently analyzed by [ML](javascript:;)-based network comparison to predict candidate stress-related genes showing expression and network differences between control and stress networks, based on 33 network topological characteristics.

## How to install mlDNA

```R
install.packages("devtools")
library(devtools)
install_github("cma2015/rsgcc")
install_github("cma2015/mlDNA")
```

## How to use mlDNA
#### Load sample data
```R
##load saltData dataset 
data(mlDNA)
##known salt stress-related genes as positive samples
positiveSamples <- as.character(sampleData$KnownSaltGenes)
##gene expression data under control condition
ControlExpMat <- as.matrix(sampleData$ControlExpMat)
##gene expression data under salt stress condition
SaltExpMat <- as.matrix(sampleData$StressExpMat)
```

#### Infer transcriptional networks from gene expression data

```R
##build transcriptional network from the first 1000 genes,
##here a higher number of cpus is suggested. 
res <- exp2net( expmat = ControlExpMat[1:1000,], method = "GCC", 
                pvalue = 0.01, cpus = 2, 
                expDescribe = "Control", connListFlag = TRUE, 
                distmatFlag = TRUE, 
                saveType = "bigmatrix", netResFileDic = netResFileDic, 
                v = rownames(ControlExpMat)[1:10],  ##for calculating distance matrix
                to =  rownames(ControlExpMat)[100:120] ) ##from "v" to "to"
```

#### Select negative samples with PSOL algorithm

```R
##generate expression feature matrix
sampleVec1 <- c(1, 2, 3, 4, 5, 6)
sampleVec2 <- c(1, 2, 3, 4, 5, 6)
featureMat <- expFeatureMatrix( expMat1 = ControlExpMat, 
                               sampleVec1 = sampleVec1, 
                               expMat2 = SaltExpMat, 
                               sampleVec2 = sampleVec2, 
                               logTransformed = TRUE, 
                               base = 2,
                               features = c("zscore", 
                                            "foldchange", "cv", 
                                            "expression"))

##positive samples
positiveSamples <- as.character(sampleData$KnownSaltGenes)
##unlabeled samples
unlabelSamples <- setdiff( rownames(featureMat), positiveSamples )

##selecting an intial set of negative samples 
##for building ML-based classification model
##suppose the PSOL results will be stored in:
PSOLResDic <- "/home/wanglab/mlDNA/PSOL/"
res <- PSOL_InitialNegativeSelection(featureMatrix = featureMat, 
                                     positives = positiveSamples, 
                                     unlabels = unlabelSamples, 
                                     negNum = length(positiveSamples), 
                                     cpus = 6, PSOLResDic = PSOLResDic)

##initial negative samples extracted from unlabelled samples with PSOL algorithm
negatives <- res$negatives

##negative sample expansion
fianlNeagatives <- PSOL_NegativeExpansion(featureMat = featureMat, positives = positiveSamples, 
                                          negatives = res$negatives, unlabels = res$unlabels, 
                                          cpus = 2, iterator = 50, cross = 5, TPR = 0.98, 
                                          method = "randomForest", plot = TRUE, trace = TRUE, 
                                          PSOLResDic = PSOLResDic,
                                          ntrees = 200 ) # parameters for ML-based classifier 
```

#### Perform cross validation

```R
##generate expression feature matrix
sampleVec1 <- c(1, 2, 3, 4, 5, 6)
sampleVec2 <- c(1, 2, 3, 4, 5, 6)
featureMat <- expFeatureMatrix( expMat1 = ControlExpMat, sampleVec1 = sampleVec1, 
                               expMat2 = SaltExpMat, sampleVec2 = sampleVec2, 
                               logTransformed = TRUE, base = 2,
                               features = c("zscore", "foldchange", "cv", "expression"))

##positive samples
positiveSamples <- as.character(sampleData$KnownSaltGenes)
##unlabeled samples
unlabelSamples <- setdiff( rownames(featureMat), positiveSamples )
idx <- sample(length(unlabelSamples))
##randomly selecting a set of unlabeled samples as negative samples
negativeSamples <- unlabelSamples[idx[1:length(positiveSamples)]]

##five-fold cross validation
seed <- randomSeed() #generate a random seed
cvRes <- cross_validation(seed = seed, method = "randomForest", 
                          featureMat = featureMat, 
                          positives = positiveSamples, 
                          negatives = negativeSamples, 
                          cross = 5, cpus = 1,
                          ntree = 100 ) ##parameters for random forest algorithm
```



## How to cite mlDNA

Chuang Ma, Mingming Xin, Kenneth A. Feldmann, Xiangfeng Wang, Machine Learning–Based Differential Network Analysis: A Study of Stress-Responsive Transcriptomes in *Arabidopsis* , *The Plant Cell*, Volume 26, Issue 2, February 2014, Pages 520–537, https://doi.org/10.1105/tpc.113.121913



## How to access help

- Comments/suggestions/bugs/issues are welcome reported [here](https://github.com/cma2015/mlDNA/issues) or contact: Chuang Ma chuangma2006@gmail.com

