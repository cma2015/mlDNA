\name{PSOL_InitialNegativeSelection}
\alias{PSOL_InitialNegativeSelection}


\title{
Initial negative set selection for building machine learning-based classification model
}


\description{
This function selects an initial negative set with the machine learning(ML)-based positive-only
 sample learning (PSOL) algorithm. The PSOL algorithm has been previously applied to predict
genomic loci encoding functional non-coding RNAs (Wang, et al. 2006). We have employed this
algorithm to identify stress-related candidate genes in Arabidopsis based on the stress 
microarray datasets (Ma and Wang, 2013).
}


\usage{
PSOL_InitialNegativeSelection(featureMatrix, positives, unlabels, 
                              negNum = length(positives), cpus = 1, PSOLResDic )
}




\arguments{
  \item{featureMatrix}{a numeric matrix recording the features for all sample.}
  \item{positives}{a character vector recording positive samples}
  \item{unlabels}{a character vector recording unlabeled samples.}
  \item{negNum}{an integer number specifying the size of negative samples will be selected.}
  \item{cpus}{an integer number specifying the number of cpus will be used for parallel computing.}
  \item{PSOLResDic}{a character string specifying the file directionry storing PSOL results.}
}



\value{
A list containing three components:
  \item{positives }{a character vector including the input positive samples.}
  \item{negatives }{a character vector recording the selected negative samples.}
  \item{unlabels }{a character vector recording the unlabeled samples.}
}



\author{
Chuang Ma and Xiangfeng Wang.
}

\references{
[1] Chunlin Wang, Chris Ding, Richard F. Meraz and Stephen R. Holbrook. PSoL: a positive sample 
only learning algorithm for finding non-coding RNA genes. Bioinformatics, 2006, 22(21): 2590-2596.

[2] Chuang Ma, Xiangfeng Wang. Machine learning-based differential network analysis: a case study
 of stress-responsive transcriptomes in Arabidopsis thaliana. 2013(Submitted).
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
                   features = c("zscore", "foldchange", 
                                 "cv","expression"))

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
                                        cpus = 6, PSOLResDic = PSOLResDic )

   ##initial negative samples extracted from unlabelled samples with PSOL algorithm
   negatives <- res$negatives

}

}



\keyword{ machine learning }
\keyword{ PSOL }
