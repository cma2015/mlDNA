\name{interactionRemoval}
\alias{interactionRemoval}



\title{
interation-removal method
}


\description{
The interation-removal method is a gene selection approach based on the removal of interesected interactions in two networks, and has been applied to analyzed the rewiring of genetic interaction maps for identying yeast genes responsive to DNA damage (Bandyopadhyay, et al., 2010; Ideker and Krogan, 2012). 
}


\usage{
interactionRemoval( adjmat1, adjmat2, threshold1, threshold2 )
}


\arguments{
  \item{adjmat1}{adjacency matrix for condition1. The adjacency matrix can be generated using adjacencymatrix function in rsgcc package.}
  \item{adjmat2}{adjacency matrix for condition2.}
  \item{threshold1}{threshold used to mask low correlations in adjmat1.}
  \item{threshold2}{threshold used to mask low correlations in adjmat2.}
}



\value{
 \item{value}{a matrix recording the specific degree of genes in the network constructed for condition 1.}
}



\author{
Chuang Ma, Xiangfeng Wang
}

\references{
[1] Sourav Bandyopadhyay, Monika Mehta, Dwight Kuo, et al. Rewiring of genetic networks in response to DNA damage.Science, 2010, 330(6009): 1385-1389.

[2] Trey Ideker and Nevan J Krogan. Differential network biology. 2012, 8:565.
}




\keyword{ network }

