\name{geneRanker}
\alias{geneRanker}


\title{
Gene Ranking
}
\description{
This function ranks genes based on differential expression analytic method.
}


\usage{
geneRanker(expmat1, expmat2, genes, 
           rankers = c("ttest", "SAM", "Limma"), verbose = FALSE)
}



\arguments{
  \item{expmat1}{a numeric matrix, gene expression matrix under condition 1.}
  \item{expmat2}{a numeric matrix, gene expression matrix under condition 2.}
  \item{genes}{a character vector, genes to be analyzed.}
  \item{rankers}{a character vector, differential expression methods}
  \item{verbose}{logical. TRUE: intermediate results will be printed to the screen.}
}


\details{
Different methods can be selected. More information can be referred in GeneSelector package 
(Boulesteix and Slawski, 2009). To run this function, please first install geneSelector package 
with the commands: source("http://bioconductor.org/biocLite.R"); biocLite("GeneSelector")
}


\value{
A list containing different components:
  \item{ranker}{a numeric matrix containing gene rank, statistic and p-value from the ranker.}
  \item{dot}{a numeric matrix containing gene rank, statistic and p-value from the rankers.}
  \item{pvalMat}{a numeric matrix containing p-values from different rankers.}
}


\references{
[1] Boulesteix A-L and Slawski M. Stability and aggregation of ranked gene lists. Brief Bioinform,
 2009, 10(5): 556-568.
}


\author{
 Chuang Ma, Xiangfeng Wang.
}



\examples{
   
\dontrun{   

   ##differential expression analysis
   res <- geneRanker(expmat1 = ControlExpMat, expmat2 = SaltExpMat, 
          genes = rownames(ControlExpMat)[1:100], 
          rankers = c("ttest", "SAM", "Limma"), 
          verbose = FALSE )
   
   ##the p-value for differential method
   res$pvalMat[1:10,]
}           
           
}


\keyword{expression}
