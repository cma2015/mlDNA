\name{optimalScore}
\alias{optimalScore}



\title{
Identifying optimal prediction score
}


\description{
The optimal prediction score is detected with the F-score measure for the machine learning-based 
classification model which can balance true positives (TPs), false positives (FPs), 
true negatives (TNs) and false negatives (FNs).
}


\usage{
optimalScore( positiveSampleScores, negativeSampleScores, beta = 2, plot = TRUE )
}

\arguments{
  \item{positiveSampleScores}{a numeric vector, the prediction scores of positive samples.}
  \item{negativeSampleScores}{a numeric vector, the prediction scores of negative samples.}
  \item{beta}{a positive numeric value, beta > 1 indicating that a higher preference is given 
  to recall than precision; beta = 1 denoting the recall and precision are weighted equally.
   beta < 1 representing that a higher preference is given to precision than recall.}
  \item{plot}{logical, TRUE indicates the distribution of F-score at different threshold of 
  prediction score is plotted. Otherwise not plotted.}
}



\value{
  A list containing two components:
  \item{statMat}{a numeric matrix recording the statistic results (i.e., prediction score, TP, FP, 
  TN, FN, Recall, TNR[TN/(TN+FP)], Precision, F-score ) at the threshold of all possible 
  prediction scores.}
  \item{optimalScore}{a numeric value, the identified optimal prediction score.}
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

    ##prediction scores of positive and negative samples from the 
    ##first round of cross validation
    positiveSampleScores <- cvRes[[1]]$positives.test.score
    negativeSampleScores <- cvRes[[1]]$negatives.test.score
    res <- optimalScore( positiveSampleScores, negativeSampleScores, 
                         beta = 2, plot = TRUE )
    
    #the optimal threshold
    res$optimalScore

    #statistic results for different threshold of prediction scores
    res$statMat[1:10,]
}
}



\keyword{ machine learning }

