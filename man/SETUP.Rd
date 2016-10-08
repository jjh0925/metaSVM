\name{SETUP}
\alias{SETUP}
\title{Setup a basic list}
\usage{SETUP(X, Y, lambda1, lambda2)
}
\arguments{
   \item{X}{The predictor variable}
   \item{Y}{The response variable}
   \item{lambda1}{A tuning parameter controls the first penalty term (group lasso)}
   \item{lambda2}{A tuning parameter controls the second penalty term (L1-lasso)}
}
\value{
\item{sm}{A list contains a basic information about model.}
}
\description{
This function generates a basic list from data sets to perform
meta-SVM.
}
\author{
SungHwan Kim, JungJun Lee, Jae-Hwan Jhong, Ja-Yong Koo
}
\seealso{\code{\link{iSVM}}}

