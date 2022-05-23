\name{ConditionSpecificGenes}
\alias{ConditionSpecificGenes}



\title{
Condition specifically expressed genes
}


\description{
This function detects condition specifcally expressed genes which are highly expressed in only
 one condition (e.g., stress), but are not or lowly expressed in the other conditions.
}


\usage{
ConditionSpecificGenes( expmat, logtransformed = TRUE, base = 2, threshold = 0.75 ) 
}



\arguments{
  \item{expmat}{a numeric matrix recording the expression level of genes at different conditions.
   The column labels are samples names. For two samples from the same condition C, their names 
   should be assigned as C.1 and C.2, respectively.}
  \item{logtransformed}{logical indicating whether the gene expression values in expmat have been
   log-transformed.}
  \item{base}{a numeric value indicating the base of logarithm.}
  \item{threshold}{a numeric value giving the threshold of condition specificity score. The 
  condition specificity score is 1, if the gene is only expressed at one condition. Otherwise,
   the condition specificity socre will be smaller than 1.}
}



\value{
	A list with following components:
    \item{CSGenes }{a character vector containing the condition specially expressed genes.}
    \item{CSScoreMat }{a data matrix recording the condition-specificity scores for all condtions,
     the maximal condition-specificity score and its corresponding condition"}
    \item{uniqueCS }{a data matrix containing the condition informaiton obtained from the
     colnames(expmat).}
}


\author{
Chuang Ma, Xiangfeng Wang
}

\references{
[1] Chuang Ma, Xiangfeng Wang. Machine learning-based differential network analysis: a study of
stress-responsive transcriptomes in Arabidopsis thaliana. 2013(Submitted).
}





\examples{

\dontrun{
  
   ##show colnames of SaltExpMat
   colnames(SaltExpMat)

   ##as the dot is used to indicate the replications in one condition, we 
   ##have to change the colname of genes at the 0.5 time point.
   colnames(SaltExpMat)[1] <- "Salt_0_5h"

   ##get condition specifically expressed genes
   res <- ConditionSpecificGenes( expmat = SaltExpMat, logtransformed = TRUE, 
                                  base = 2, threshold = 0.75 ) 

   ##number of condition specifically expressed genes
   length(res$CSGenes)
 }  

}



\keyword{ statistic }

