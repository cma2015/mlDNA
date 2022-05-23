\name{expFeatureMatrix}
\alias{expFeatureMatrix}



\title{
Expression-based feature matrix
}
\description{
This function generates expression-based features for each gene with the consideration of z-scores,
 fold changes and actural expression values.
}



\usage{
expFeatureMatrix(expMat1, sampleVec1, expMat2, sampleVec2, 
                 logTransformed = TRUE, base = 2, 
                 features = c("zscore", "foldchange", "cv", "expression"))
}




\arguments{
  \item{expMat1}{a numeric matrix recording gene expression data from condition 1. Each row 
  represent the expression values of one gene, each column represents the expression  values of all
   genes from one biological experiment.}
  \item{sampleVec1}{a numeric vector representing biological replication (or different time points)
   and technical replication for condition 1. For instance, c(1,1,2,2,3,3) denotes three biological
    replications, and two technicial replications for each biological sample.  }
  \item{expMat2}{a numeric matrix recording gene expression data from condition 2.}
  \item{sampleVec2}{a numeric vector representing biological experiments and technical replications
   for condition 2. Definition is similar as sampleVec1. }
  \item{logTransformed}{logic value, TRUE indicates that the expression data in expMat1 and expMat2
   have been log-transformed. }
  \item{base}{base of log transformation. }
  \item{features}{A character fector specifying different features will be used to generate
   feature matrix. }
}



\value{
  \item{value}{a numeric matrix with express-based features for each gene.}
}


\author{
   Chuang Ma, Xiangfeng Wang
}



\note{
The expression-based features including:

(1) expression variance across all samples calculated with z-scores of genes under two conditions;

(2) fold change of gene expression under two conditions;

(3) expression values of genes in two conditions used for retaining the information from 
actual expression values. 
}



\examples{

 \dontrun{
 
   ##generated expression features
   sampleVec1 <- c(1, 2, 3, 4, 5, 6)
   sampleVec2 <- c(1, 2, 3, 4, 5, 6)
   featureMat <- expFeatureMatrix( expMat1 = ControlExpMat, sampleVec1 = sampleVec1, 
                                   expMat2 = SaltExpMat, sampleVec2 = sampleVec2, 
                                   logTransformed = TRUE, base = 2, 
                               features = c("zscore", "foldchange", "cv", "expression"))
}

}

\keyword{differential expression}

