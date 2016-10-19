\name{iSVM}
\alias{iSVM}
\title{Integrated Support Vector Machine}
\usage{iSVM(sm, is.constant=TRUE)}
\arguments{
   \item{sm}{A list generated from \code{SETUP}}
   \item{is.constant}{if TRUE, an intercept is included in the model; default is TRUE.}
}
\value{A list contains information on the final model}
\description{
This function develops an estiamation and variable selection algorithm for
meta-analytic framework based on the support vector machine (Meta-SVM).
}
\details{
The proposed Meta-SVM is motivated by the recent meta-analytic method exploiting the logistic regression.
We develop a novel implementation strategy in spirit of Newton’s method in the Meta-SVM. For the most part, the objective function of SVM is formed with the hinge loss and a range of penalty terms (e.g., L1-lasso, group lasso and etcs).
We particularly adopts the sparse group lasso enabling to capture both common and study specific genetic signals among all studies.
}
\author{
SungHwan Kim, JungJun Lee, Jae-Hwan Jhong, Ja-Yong Koo
}
\seealso{\code{\link{SETUP}}}
\examples{
library(metaSVM)
fpath = system.file("extdata", "DList_lung_svm.Rdata", package="metaSVM")
DList = get(load(fpath))
Y = list()
for (k in 1:length(DList))
{
   tmp = as.numeric(rownames(DList[[k]]) == "1")
   tmp[tmp == 0] = -1
   Y[[k]] = tmp
}
for(i in 1:length(DList))
   rownames(DList[[i]]) = 1:dim(DList[[i]])[1]
sm = SETUP(X = DList, Y = Y, lambda1 = 0.1, lambda2 = 0.1)
Res = iSVM(sm, is.constant = TRUE)
}



