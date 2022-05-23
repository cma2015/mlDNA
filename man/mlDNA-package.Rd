\name{mlDNA-package}
\alias{mlDNA-package}
\alias{mlDNA}
\docType{package}

\title{
Machine Learning-based Differential Network Analysis
}
\description{
This package provides functions for performing the machine learning (ML)-based differential 
network analysis of transcriptome data. It can be used to: 

1) perform machine learning-based gene filtering with positive sample-only learning algorithm for 
identifying a set of candidate genes with four different classes of expression characteristics, 
including the absoulte expression values, the within-condition expression variations, 
the between-condition expression changes, and the coefficient of variation;

2) construct gene co-expression networks from gene expression data with seven optional methods,
including five correlation and two non-correlation measures;

3) perform a comprehensive network comparision with more than thirty network-based characteristics
 including degree, closeness, eccentricity, eigenvector, and PageRank;

4) identify biologically important genes with different ML algorithms by combining network-based 
characteristics generated from differential network analysis; 

5) estimate the covergence degree between different experimential conditions;

7) quantify the activity of pathways;

8) detect condition specifcally expressed genes.

The tutorial of the mlDNA package can be found at: http://www.cmbb.arizona.edu/mlDNA.
}


\details{
\tabular{ll}{
Package: \tab mlDNA\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-11-18\cr
License: \tab GPL(>=2)\cr
}

}


\author{
Chuang Ma, Xiangfeng Wang

Maintainer: Chuang Ma <chuangma2006@gmail.com>
}


\references{
[1] Chuang Ma, Xiangfeng Wang. Machine learning-based differential network analysis: a study of
stress-responsive transcriptomes in Arabidopsis thaliana. 2013(Submitted).
}



\keyword{ package }

