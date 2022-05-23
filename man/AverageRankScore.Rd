\name{AverageRankScore}
\alias{AverageRankScore}



\title{
Average-based rank score
}


\description{
This function calculates the activities of pathways in the whole genome with the average-based
 rank scoring algorithm (Yang, et al., 2011; Ma and Wang, 2013). This rank-based statistics is 
 robust for directly comparing the activities of pathways with different gene numbers under 
 different experimental conditions, since it produces a normalized value with the consideration
  of gene number in the analyzed pathways and whole genomes.
}


\usage{
AverageRankScore( featureMat, selGenes )
}


\arguments{
  \item{featureMat}{a numeric matrix recording the expression levels or changes of all genes in
   the genome at given conditions.}
  \item{selGenes}{a character vector recording a set of genes in the analyzed pathway.}
}



\value{
 \item{value}{a numeric vector recording the activities of interested genes (selGenes)
  at different conditions.}
}



\author{
Chuang Ma, Xiangfeng Wang
}

\references{
[1] Huang Yang, Chao Cheng and Wei Zhang. Average rank-based score to measure deregulation
 of molecular pathway gene sets. PLoS One, 2011, 6(11): e27579.

[2] Chuang Ma, Xiangfeng Wang. Machine learning-based differential network analysis: a study of
stress-responsive transcriptomes in Arabidopsis thaliana. 2013(Submitted).
}


\examples{
  
\dontrun{

   ##generate expression feature matrix
   sampleVec1 <- c(1, 2, 3, 4, 5, 6)
   sampleVec2 <- c(1, 2, 3, 4, 5, 6)
   featureMat <- expFeatureMatrix( expMat1 = ControlExpMat, sampleVec1 = sampleVec1, 
                                   expMat2 = SaltExpMat, sampleVec2 = sampleVec2, 
                                   logTransformed = TRUE, base = 2,
                                   features =  "foldchange" )

   ##for an interested set of genes, the average-based rank score can be calculated:
   genes <- rownames(featureMat)[1:100]
   res <- AverageRankScore( featureMat = featureMat, selGenes = genes )

}

}



\keyword{ statistic }

