\name{classifier}
\alias{classifier}



\title{
Building machine learning-based classification models
}


\description{
This function builds classification models with different machine learning algorithms including
 random forest (randomForest), support vector machine (svm), and neural network (nnet).
}
\usage{
classifier( method = c("randomForest", "svm", "nnet" ), 
            featureMat, positiveSamples, negativeSamples, 
            tunecontrol = tune.control(sampling = "cross", cross = 5), ...)
}



\arguments{
  \item{method}{a character string specifying machine learning method used to construct 
  prediction models. Currently the user-optional values are "randomForest", "svm" or "nnet".}
  \item{featureMat}{a numeric matrix recording feature values for each sample.}
  \item{positiveSamples}{a character vector recording positive samples used to
   train classification model.}
  \item{negativeSamples}{a character vector recording negative samples used to
   train classification model. }
  \item{tunecontrol}{control parameters for the tune function. The value is assigned by the
   function "tune.control" in e1071 package. }
  \item{\dots}{Further parameters passed to train classification model. }
}


\details{

  For the random forest algorithm, the important parameters are mtry (number of features randomly
   selected for bulding the decision tree. Default:sqrt(ncol(featureMat))) and ntree (number of 
   trees to be built. Default:500). More information about the parameters related to random forest
    can be found in R package RandomForest.

  For the svm algorithm, the default kernel function is the "radial"(radial basis function; RBF).
   Other optional kernels are the "linear", the "polynomial", and the "sigmoid" fuctions. The range
    of degree (for the kernel type of polynomial), gamma (for radial), coef0 (for polynomial 
    and sigmoid) and cost (the cost parameters for all kernels) should be designated. Please refer
   the description of "svm" function in R package e1071 for more information 
   about the parameters of svm. 

  For the neural network, the important parameters are size (number of units in the hidden layer),
   decay (parameter for weight decay. Default 0), trace (switch for tracing optimization. 
   Default TRUE.). More information about the parameters of neural network is described in
    nnet and e1071 packages.

}


\value{
An object of class tune (see tune function in e1071 package) is returned, including the components:

 \item{best.parameters }{a 1 x k data frame, k number of parameters.}
 \item{best.performance }{best achieved performance.}
 \item{best.model }{ the model trained on the complete training data 
 using the best parameter combination.}

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
   ##randomly selecting 1000 unlabeled samples as negative samples
   negativeSamples <- unlabelSamples[idx[1:1000]]

   ##for random forest, and using five-fold cross validation for obtaining optimal parameters
   cl <- classifier( method = "randomForest", 
                     featureMat = featureMat, 
                     positiveSamples = positiveSamples, 
                     negativeSamples = negativeSamples,
                     tunecontrol = tune.control(sampling = "cross", cross = 5),
                     ntree = 100 ) #build 100 trees for the forest


   ##svm and using five-fold cross validation for obtaining optimal parameters
   cl <- classifier( method = "svm", featureMat = featureMat, 
                     positiveSamples = positiveSamples, 
                     negativeSamples = negativeSamples,
                     tunecontrol = tune.control(sampling = "cross", cross = 5),
                     kernel = "radial", 
                     probability = TRUE,
                     ranges = list(gamma = 2^(-2:2), 
                      cost = 2^(-4:4)) ) #for radial kernel and the parameter space.

   ##neural network, using one split for training/validation set
   cl <- classifier( method = "nnet", featureMat = featureMat, 
                     positiveSamples = positiveSamples, 
                     negativeSamples = negativeSamples,
                     tunecontrol = tune.control(sampling = "fix"), 
                     trace = TRUE, size = 10 ) #for nnet parameters.
  

}
}


\keyword{machine learning}
