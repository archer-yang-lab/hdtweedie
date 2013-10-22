\name{coef.cv.HDtweedie}
\alias{coef.cv.HDtweedie}
\title{get coefficients or make coefficient predictions from a "cv.HDtweedie" object.}
\description{
This function gets coefficients or makes coefficient predictions from a cross-validated \code{HDtweedie} model,
using the \code{"cv.HDtweedie"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{coef}{cv.HDtweedie}(object,s=c("lambda.1se","lambda.min"),...)
}
\arguments{
	\item{object}{fitted \code{\link{cv.HDtweedie}} object.}
	\item{s}{value(s) of the penalty parameter \code{lambda} at which
	predictions are required. Default is the value \code{s="lambda.1se"} stored
	on the CV \code{object}, it is the largest value of \code{lambda} such that error is
	within 1 standard error of the minimum. Alternatively \code{s="lambda.min"} can be
	used, it is the optimal value of \code{lambda} that gives minimum
	cross validation error \code{cvm}. If \code{s} is numeric, it is taken as the value(s) of
	\code{lambda} to be used.}
	\item{\dots}{not used. Other arguments to predict.} 
}
\details{This function makes it easier to use the results of
cross-validation to get coefficients or make coefficient predictions.}
\value{The coefficients at the requested values for \code{lambda}.}
\author{Wei Qian, Yi Yang and Hui Zou\cr
Maintainer: Wei Qian  <weiqian@stat.umn.edu>}
\references{
Qian, W., Yang, Y., Yang, Y. and Zou, H. (2013), ``Tweedie's Compound
    Poisson Model With Grouped Elastic Net,'' submitted to \emph{Journal of Computational and Graphical Statistics}.\cr

Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.HDtweedie}}, and \code{\link{predict.cv.HDtweedie}} methods.}
\examples{
# load HDtweedie library
library(HDtweedie)

# load data set
data(auto)

# 5-fold cross validation using the lasso
cv0 <- cv.HDtweedie(x=auto$x,y=auto$y,p=1.5,nfolds=5)

# the coefficients at lambda = lambda.1se
coef(cv0)

# define group index
group1 <- c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)

# 5-fold cross validation using the grouped lasso 
cv1 <- cv.HDtweedie(x=auto$x,y=auto$y,group=group1,p=1.5,nfolds=5)

# the coefficients at lambda = lambda.min
coef(cv1, s = cv1$lambda.min)
}
\keyword{models}
\keyword{regression}
