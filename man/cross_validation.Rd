\name{cross_validation}
\alias{cross_validation}


\title{
Cross validation method
}


\description{
The ML-based classification model is trained and tested with N-fold cross validation method.
}


\usage{
cross_validation(seed = 1, method = c("randomForest", "svm", "nnet" ), 
                 featureMat, positives, negatives, cross = 5, cpus = 1, ...)
}




\arguments{
  \item{seed}{an integer number specifying a random seed for randomly partitioning dataset. }
  \item{method}{a character string specifying machine learning method. Possible values 
  are "randomForest", "nnet" or "svm"}
  \item{featureMat}{a numeric feature matrix.}
  \item{positives}{a character vector reocrding positive samples}
  \item{negatives}{a character vector recording negative samples.}
  \item{cross}{number of fold for cross validation.}
  \item{cpus}{an integer number specifying the number of cpus to be used for parallel computing.}
  \item{\dots}{Further parameters used to cross validation. Same with the parameters 
  used in the classifer function.}
}


\details{
In machine learning, the cross validation method has been widely used to evaluate the performance
 of ML-based classification models (classifiers). 

For N-fold cross validation, positive and negative samples are randomly partitioned into N groups 
with approximately equal amount of samples, and each group is successively used for testing the 
performance of the ML-based classifier trained with the other N-1 groups of positive and negative
 samples. 

For each round of cross validation, the prediction accuracy of the ML-based classifier was assessed
 using the receiver operating characteristic (ROC) curve analysis.The ROC curve is a two-dimensional
  plot of the false positive rate (FPR, x-axis) against the true positive rate (TPR, y-axis) at all
   possible thresholds. The value of area under the ROC curve (AUC) was used to quantitatively score
    the prediction accuracy of the ML-based classifer. The AUC value is ranged from 0 to 1.0, with 
    higher AUC value indicates a better prediction accuracy of the ML-based classifer. 

After N groups have been successively used as the testing set, the N sets of (FPR, TPR) pairs were 
imported into R package ROCR to visualize the ROC curves. The mean value of N AUCs was then computed
 as the overall performance of the ML-based classification model.
}


\value{
A list recording results from each fold cross validation including the components:
  \item{positves.train }{positive samples used to train prediction model. }
  \item{negatives.train }{negative samples used to train prediction model. }
  \item{positives.test }{positive samples used to test prediction model. }
  \item{negatives.test }{negative samples used to test prediction model. }
  \item{ml }{machine learning method. }
  \item{classifier }{prediction model constructed with the best parameters obtained 
  from training dataset. }
  \item{positives.train.score }{scores of postive samples in training dataset 
  predicted by classifier. }
  \item{positives.train.score }{scores of postive samples in training dataset 
  predicted by classifier. }
  \item{positives.test.score }{scores of postive samples in testing dataset 
  predicted by classifier. }  
  \item{negatives.test.score }{scores of negative samples in testing dataset 
  predicted by classifier. }  
  \item{train.AUC }{AUC value of the ML-based classifer on training dataset. }  
  \item{test.AUC }{AUC value of the ML-based classifer on testing dataset. } 
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
                             positives = positiveSamples, 
                             negatives = negativeSamples, 
                             cross = 5, cpus = 1,
                             ntree = 100 ) ##parameters for random forest algorithm

   ##get AUC values for five rounds of cross validation
   aucVec <- rep(0, 5) 
   for( i in 1:5 ) 
     aucVec[i] = cvRes[[i]]$test.AUC
  
   
   ##average AUC values as the final performance of the ML-based classifier
   mean(aucVec)

 
}

}

\keyword{ machine learning }

