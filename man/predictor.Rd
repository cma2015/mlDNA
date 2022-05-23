\name{predictor}
\alias{predictor}


\title{
Prediction
}


\description{
Score each sample with the machine learning-based classification prediction model alreadly trained 
with training dataset.
}


\usage{
predictor(method = c("randomForest", "svm", "nnet" ), classifier, featureMat)
}



\arguments{
  \item{method}{character string specifying the machine learning algorithm used to 
  buld classification model.}
  \item{classifier}{trained prediction model obtained from the classifier function.}
  \item{featureMat}{a numeric matrix; feature matrix containing samples to be scored 
  and their feature values. }
}


\value{
  \item{value }{A numeric vector containing the prediction score of each sample.}
}


\author{
Chuang Ma, Xiangfeng Wang
}


\examples{

\dontrun{

   ##generate expression feature matrix
   sampleVec1 <- c(1, 2, 3, 4, 5, 6)
   sampleVec2 <- c(1, 2, 3, 4, 5, 6)
   featureMat <- expFeatureMatrix( 
           expMat1 = ControlExpMat, sampleVec1 = sampleVec1, 
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

   ##for random forest, and using five-fold cross validation 
   ##for obtaining optimal parameters
   cl <- classifier( method = "randomForest", featureMat = featureMat, 
        positiveSamples = positiveSamples, negativeSamples = negativeSamples,
        tunecontrol = tune.control(sampling = "cross", cross = 5), 
        ntree = 100 ) #build 100 trees for the forest


   ##constructed prediction model
   predModel <- cl$best.model

   ##perform prediction
   predResult <- predictor(method = "randomForest", 
                           classifier = predModel, 
                           featureMat = featureMat)

}

}


\keyword{machine learning}
