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
\item{sm}{A list contains all related objects to implement the Meta-SVM.}
}
\description{
This function generates a list containing
all relevant objects to perform the meta-SVM.
}
\author{
SungHwan Kim, Jae-Hwan Jhong, JungJun Lee, Ja-Yong Koo
}
\seealso{\code{\link{iSVM}}}

