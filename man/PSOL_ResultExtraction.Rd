\name{PSOL_ResultExtraction}
\alias{PSOL_ResultExtraction}


\title{
PSOL result extraction
}

\description{
This function extracts the PSOL result.
}


\usage{
PSOL_ResultExtraction(PSOLResDic, iterations = c(1:4) )
}



\arguments{
  \item{PSOLResDic}{PSOL result file directory. }
  \item{iterations}{a numeric vector specifying the results at given iteration times 
  will be extracted. }
}



\value{
A list with length(iterations) components. For each components, 
there is a list with four components: 
  \item{AUC}{The AUC from cross validation experiments testing classifiers 
  with positives and negatives.}
  \item{positives}{positive samples}
  \item{negatives}{negative samples}
  \item{unlabels}{unlabeled samples}
}

\author{
Chuang Ma, Xiangfeng Wang
}



\examples{
   
\dontrun{

   ##after run PSOL_InitialNegativeSelection and PSOL_NegativeExpansion
   ##extract the PSOL results at specified iteration times with the command:
   PSOLResDic <- "/home/wanglab/mlDNA/PSOL/"
   PSOLRes <- PSOL_ResultExtraction( PSOLResDic = PSOLResDic, iterations = c(1:4) )
}
}


\keyword{ machine learning }
\keyword{ PSOL }
