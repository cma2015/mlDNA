\name{exp2net}
\alias{exp2net}


\title{
Inferring transcriptional networks from gene expression data
}


\description{
This function infers transcriptional networks from gene expression data with different 
statistical methods, including five correlation measures (i.e., the Gini correlation 
coefficient [GCC], the Pearson's product moment correlation coefficient [PCC], the  Kendall tau 
rank correlation coefficient [KCC], the Spearman's rank correlation coefficient [SCC] and the 
Tukey's biweight correlation coefficient [BiWt]) and two non-correlation measures (mutual 
information [MI] and the maximal information-based nonparametric exploration [MINE]).
}


\usage{
exp2net( expmat, method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE"),
         pvalue = 0.01, cpus = 1, expDescribe = "Control", 
         connListFlag = TRUE, distmatFlag = TRUE, saveType = "bigmatrix",
         netResFileDic, ...  )
}



\arguments{
  \item{expmat}{a numberic matrix recording gene expression data.}
  \item{method}{a character string specifying the statistical method will be used to calculating 
  the associations between any pairs of genes.}
  \item{pvalue}{a numeric value denoting the significance level of the association will be used to
   filter unsignficant interactions (i.e., edge) in the network.}
  \item{cpus}{an integer specifying the number of cpus will be used for parallel computing.}
  \item{expDescribe}{an character string describing the expmat.}
  \item{connListFlag}{a logical value indicating whether the connected genes for each gene 
  will be recorded.}
  \item{distmatFlag}{a logical value indicating whether the distance matrix will be calculated.}
  \item{saveType}{an character string indicating the format ("matrix", "bigmatrix") of matrix.}
  \item{netResFileDic}{a character string specifying the file directory will be used to store
   network-related results.}
  \item{\dots}{Furture parameters for calcluating distances between two gene sets. For instance,
   v = c(g1, g2, ..., gn), to = c(g1, g3, ..., gm). }
}


\value{
A list with 12 components: 
  \item{expmat}{the input gene expression data.}
  \item{method}{the method used to calcluate the association between two genes.}
  \item{pvalue}{the significance level used to detect edges in the network.}
  \item{expDescribe}{the characterized string for gene expression data.} 
  \item{netResFileDic}{the file directory for storing network-related result.  } 
  \item{adjmat}{adjacency matrix recording the association between any pairs of genes
   in the big.matrix format.} 
  \item{adjmat_backingfile}{the root name for the file for the cache of adjmat.
   Default: expDescribe_method_adjmat_bfile} 
  \item{adjmat_descriptorfile}{the file to be used for the description of the adjmat.
   Default: expDescribe_method_adjmat_dfile}
  \item{threshold}{the correlation score at the significance level of pvalue.} 
  \item{graph}{an igraph object for the constructed network in the edgelist format.
   This object is save in the file: expDescribe_graph.} 
  \item{connectivityList}{a list; For each component, it is a list with three component:
   "pos" (connected genes with positive correlations), 
   "neg" (connected genes with negative correlations), "all" (all connected genes)} 
  \item{distmatrix}{a numeric matrix; the shorest-path distance between any pair of 
  genes in the network.} 
 
}

\author{
Chuang Ma, Xiangfeng Wang.
}


\note{
[1] The GCC, PCC, SCC and KCC calcluate the adjacency matrix more quickly than BiWt, MI and MINE.

[2] The threshold is determined with the permutation method by generating the background 
distribution of correlations by permuting the expression levels of nrow(expmat) genes from the
 original expression dataset(expmat) (Carter et al., 2004).

[3] The adjacency and distance matrix can be stored in big.matrix format which can be used to 
greatly save the memory space. However, this big.matrix optional can only be used on Linux. 
More information about the big.matrix can be found in the R package bigmemory.

[4] The functions for graph analysis (i.e., getting information of nodes and edges) can be found
 in the R package igraph.

[5] The calculation of distance matrix is time-consuming for large-scale network.
}

\references{
[1] Scott L. Carter, Christian M. Brechbuhler, Michael Griffin and Andrew T. Bond. Gene 
co-expression network topology provides a framework for molecular characterization of cellular
 state. Bioinformatics, 2004, 20(14): 2242-2250.
}




\examples{

 \dontrun{

   ##suppose the network-related results are stored at:
   netResFileDic = "/home/wanglab/mlDNA/network/"
   
   ##build transcriptional network from the first 1000 genes,
   ##here a higher number of cpus is suggested. 
   res <- exp2net( expmat = ControlExpMat[1:1000,], method = "GCC", 
                   pvalue = 0.01, cpus = 2, 
                   expDescribe = "Control", connListFlag = TRUE, 
                   distmatFlag = TRUE, 
                   saveType = "bigmatrix", netResFileDic = netResFileDic, 
                   v = rownames(ControlExpMat)[1:10],  ##for calculating distance matrix
                   to =  rownames(ControlExpMat)[100:120] ) ##from "v" to "to"
                   
}
}



\keyword{network}

