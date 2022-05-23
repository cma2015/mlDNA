\name{PSOL_NegativeExpansion}
\alias{PSOL_NegativeExpansion}


\title{
 Negative expansion
}
\description{
This function expands the negative sample set using PSOL algorithm.
}


\usage{
PSOL_NegativeExpansion(featureMat, positives, negatives, unlabels, cpus = 1, 
                       iterator = 50, cross = 5, TPR = 0.98, method = "randomForest", 
                       plot = TRUE, trace = TRUE, PSOLResDic, ...)
}




\arguments{
  \item{featureMat}{a feature matrix recording the feature values for all samples.}
  \item{positives}{a character string recording the positive samples.}
  \item{negatives}{a character string recording the negative samples.}
  \item{unlabels}{a character string recording the unlabeled samples.}
  \item{cpus}{an integer value, cpu number}
  \item{iterator}{an integer value, iterator times.}
  \item{cross}{an integer value, cross-times cross validation.}
  \item{TPR}{a numeric value ranged from 0 to 1.0, used to select the prediction score cutoff.}
  \item{method}{a character string, machine learing method}
  \item{plot}{a logic value specifies whether the score distribution of positive and 
  unlabeled samples will be plotted.}
  \item{trace}{logic. TRUE: the intermediate results will be saved as ".RData" format.}
  \item{PSOLResDic}{a character string, PSOL Result directory}
  \item{\dots}{Further parameters used in PSOL_ExpandSelection. 
  see the further parameters in function classifier.}
}



\value{
The PSOL-related results are output in the "resultDic" directory.
}


\author{
 Chuang Ma, Xiangfeng Wang.
}


\examples{

\dontrun{

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
   PSOL_NegativeExpansion(featureMat = featureMat, positives = positiveSamples, 
                          negatives = res$negatives, unlabels = res$unlabels, 
                          cpus = 2, iterator = 50, cross = 5, TPR = 0.98, 
                          method = "randomForest", plot = TRUE, trace = TRUE, 
                          PSOLResDic = PSOLResDic,
                          ntrees = 200 ) # parameters for ML-based classifier 

}

}


\keyword{machine learning}
\keyword{ PSOL }

