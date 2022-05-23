\name{plotROC}
\alias{plotROC}



\title{
Plotting receiver operating characteristic(ROC) Curves
}


\description{
This function plots ROC curves for estimating the performance of machine learning-based 
classification model in cross validation experiments.
}


\usage{
plotROC(cvRes)
}

\arguments{
  \item{cvRes}{results from the "cross_validation" function.}
}



\value{
 A ROC plot
}


\author{
Chuang Ma, Xiangfeng Wang
}


\examples{
 \dontrun{

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
                             positives = positiveSamples, negatives = negativeSamples, 
                             cross = 5, cpus = 1,
                             ntree = 100 ) ##parameters for random forest algorithm

   #plot ROC curve
   plotROC(cvRes)
}
}



\keyword{ machine learning }

