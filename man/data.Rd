\name{data}
\alias{data}
\alias{saltData}
\non_function{}


\title{
	example data for mlDNA 
}
\description{
	salt stress expression data and 50 known salt-related genes
}
\usage{
	data(mlDNA)
}


\examples{

   ##load saltData dataset 
   data(mlDNA)

   ##known salt stress-related genes as positive samples
   positiveSamples <- as.character(sampleData$KnownSaltGenes)

   ##gene expression data under control condition
   ControlExpMat <- as.matrix(sampleData$ControlExpMat)

   ##gene expression data under salt stress condition
   SaltExpMat <- as.matrix(sampleData$StressExpMat)
   
}


\keyword{
 datasets
}
