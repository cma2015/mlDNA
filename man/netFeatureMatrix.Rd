\name{netFeatureMatrix}
\alias{netFeatureMatrix}


\title{
Generating network feature matrix
}


\description{
This functions generates a feature matrix containing 33 network characteristics by 
differential network analysis.
}


\usage{

netFeatureMatrix( net1, net2, nodes = NULL, knodes = NULL, 
                  cpus = 1, verbose = TRUE, netResFileDic,
                  features = c( "expDistance", "ASC", "corDistance", 
                                "AllConnectivity", "PosConnectivity", 
                                 "NegConnectivity", "closeness", 
                                 "eccentricity", "eigenvector", 
                                 "page.rank", "dis2knodes", 
                                 "closeness2knodes", "eccenticity2knodes") )
}



\arguments{

  \item{net1}{exp2net output for condition 1.}
  \item{net2}{exp2net output for condition 2.}
  \item{nodes}{a character vector, a set of genes will be analyzed.}
  \item{knodes}{a character vector, a set of interested genes (e.g., known stress-related genes).}
  \item{cpus}{an integrator value, the number of cpus will be used for parallel computing.}
  \item{verbose}{logical value indicate whether the progress information will be output.}
  \item{netResFileDic}{file directory where the network-related results will be output.}
  \item{features}{a character vector specifying the network characteristics will be analyzed.}
}


\note{
1) More information about these network properties has been described in (Ma and Wang, 2013).

2) To run netFeatureMatrix, the parameters "v" and "to" in "exp2net" function should be the 
whole genes in expmat.
}


\value{
a numeric matrix, feature matrix generated from network comparision analysis.
}



\author{
Chuang Ma, Xiangfeng Wang.
}


\references{
[1] Chuang Ma, Xiangfeng Wang. Machine learning-based differential network analysis:a case 
study of stress-responsive transcriptomes in Arabidopsis thaliana. 2013(Submitted).
}





\examples{



\dontrun{
   
   ##suppose the network-related results are stored at:
   netResFileDic = "/home/wanglab/mlDNA/network/"
   
   ##only consider a subset of genes here
   genes <- unique ( c( rownames(ControlExpMat)[1:1000], positiveSamples[1:100] ) )
   res_ControlSub <- exp2net( expmat = ControlExpMat[genes,], method = "GCC", 
                              pvalue = 0.01, cpus = 6, expDescribe = "Control_Sub", 
                              connListFlag = TRUE, distmatFlag = TRUE, 
                              saveType = "bigmatrix", netResFileDic = netResFileDic, 
                              v = genes, to = genes )
                   
    res_StressSub <- exp2net( expmat = SaltExpMat[genes,], method = "GCC", 
                              pvalue = 0.01, cpus = 6, expDescribe = "Stress_Sub", 
                              connListFlag = TRUE, distmatFlag = TRUE, 
                              saveType = "bigmatrix", netResFileDic = netResFileDic, 
                              v = genes, to = genes )

    ##generate network feature matrix
    nodes <- genes
    knodes <- intersect( genes, positiveSamples )
    netFeatureMat <- netFeatureMatrix( net1 = res_ControlSub, net2 = res_StressSub, 
                                       nodes = nodes, knodes = knodes, 
                                       cpus = 2, verbose = TRUE, 
                                       netResFileDic = netResFileDic,
                                       features = c( "expDistance", "ASC", 
                                             "corDistance", "AllConnectivity", 
                                             "PosConnectivity", "NegConnectivity", 
                                             "closeness","eccentricity", 
                                             "eigenvector", "page.rank", 
                                             "dis2knodes", "closeness2knodes", 
                                             "eccenticity2knodes") )
                                                     
}


}


\keyword{network}

