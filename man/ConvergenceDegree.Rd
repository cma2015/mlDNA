\name{ConvergenceDegree}
\alias{ConvergenceDegree}



\title{
Calculating covergence degree
}


\description{
This function calculates the covergence degree between two sets of elements, and has been used to 
estimate the degree of covergence between genes responsive to different stresses in Arabidopsis
 (Ma and Wang, 2013)
}


\usage{
ConvergenceDegree( vec1, vec2 )
}

\arguments{
  \item{vec1}{a vector recording a set of elements for one condition.}
  \item{vec2}{a vector recording a set of elements for the other condition.}
}



\value{
 \item{value}{an numeric value}
}


\author{
Chuang Ma, Xiangfeng Wang
}

\references{
[1] Chuang Ma, Xiangfeng Wang. Machine learning-based differential network analysis: a study of
stress-responsive transcriptomes in Arabidopsis thaliana. 2013(Submitted).
}




\examples{
   vec1 <- c(1, 2, 4, 10)
   vec2 <- c(1, 10, 15:20)
   ConvergenceDegree( vec1, vec2 )
}



\keyword{ statistic }

