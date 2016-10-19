\name{iSVM}
\alias{iSVM}
\title{Integrated Support Vector Machine}
\usage{iSVM(sm, is.constant=TRUE)}
\arguments{
   \item{sm}{A list generated from \code{SETUP}}
   \item{is.constant}{if TRUE, an intercept is included in the model.}
}
\value{A list contains information on the final model}
\description{
This function develops an estimation and variable selection algorithm for
meta-analytic framework based on the meta-analytic support vector machine (Meta-SVM).
}
\details{
The proposed Meta-SVM is motivated by the recent meta-analytic method exploiting the logistic regression.
We develop a novel implementation strategy in spirit of Newton’s method in the Meta-SVM. For the most part, the objective function of SVM is formed with the hinge loss and a range of penalty terms (e.g., L1-lasso, group lasso and etcs).
We particularly adopts the sparse group lasso enabling to capture both common and study specific genetic signals among all studies.
}
\author{
SungHwan Kim, Jae-Hwan Jhong, JungJun Lee, Ja-Yong Koo
}
\seealso{\code{\link{SETUP}}}
\examples{
library(metaSVM)
fpath = system.file("extdata", "DList_lung_svm.Rdata", package = "metaSVM")
DList = get(load(fpath))

#==============#
# The following stage proceeds to rename the binary outcome variable (1:IPF/0:control)
# and gene features.
#==============#
Y = list()
for (k in 1:length(DList))
{
   tmp = as.numeric(rownames(DList[[k]]) == "1")
   tmp[tmp == 0] = -1
   Y[[k]] = tmp
}
for(i in 1:length(DList))
   rownames(DList[[i]]) = 1:dim(DList[[i]])[1]

#==============#
# As a part of preprocessing, SETUP undertakes tasks to collect all required data
# and to put them into the list form.
#==============#
sm = SETUP(X = DList, Y = Y, lambda1 = 0.1, lambda2 = 0.1)
# Fit the meta-SVM model using the example data set
# The following process takes about 25 minutes on average.
# Res = iSVM(sm, is.constant = TRUE)
}
\format{
Dlist_lung_svm.Rdata consists of three microarray profiles
of the idiopathic pulmonary fibrosis (IPF) gathered from
KangA (Kang et al., 2012), Konishi (Konishi et al., 2009),
and Pardo (Pardo et al., 2005),
each containing 63, 38, and 24 samples, respectively.
For pre-processing, we performed unbiased filtering via
standard deviation such that it remains 39 features.
}
\source{
The time required to fit the metaSVM model was measured on R
version 3.3.1 (2016-06-21) running on macOS Sierra,
in a laptop with an Intel Core i7 CPU (2.60 GHz) and 16.0 GB RAM.
}
\references{
Kang, D., Sibille, E., Kaminski, N., Tseng, G. MetaQC: objective quality control and inclusion/exclusion criteria for genomic meta-analysis. Nucleic Acids Res., 40, e15, 2012.

Konishi K., Gibson K., Lindell K., Richards T., Zhang Y., et al. Gene expression profiles of acute exacerbations
of idiopathic pulmonary fibrosis Am J Respir Crit Care Med, 180(2):167–75, 2009.

Pardo A., Selman M. Role of matrix metaloproteases in idiopathic pulmonary fibrosis Fibrogenesis Tissue Repair, 5(Suppl 1):S9, 2012.
}

