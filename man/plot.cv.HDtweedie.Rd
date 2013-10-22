\name{plot.cv.HDtweedie}
\alias{plot.cv.HDtweedie}
\title{plot the cross-validation curve produced by cv.HDtweedie}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv} function from the \code{glmnet} package.}
\usage{
\method{plot}{cv.HDtweedie}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{fitted \code{\link{cv.HDtweedie}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or
	its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}
\details{A plot is produced.}

\author{Wei Qian, Yi Yang and Hui Zou\cr
Maintainer: Wei Qian  <weiqian@stat.umn.edu>}
\references{
Qian, W., Yang, Y., Yang, Y. and Zou, H. (2013), ``Tweedie's Compound
    Poisson Model With Grouped Elastic Net,'' submitted to \emph{Journal of Computational and Graphical Statistics}.\cr

Friedman, J., Hastie, T., and Tibshirani, R. (2010), ``Regularization paths for generalized
linear models via coordinate descent,'' \emph{Journal of Statistical Software}, 33, 1.\cr
\url{http://www.jstatsoft.org/v33/i01/}
}
\seealso{\code{\link{cv.HDtweedie}}.}
\examples{
# load HDtweedie library
library(HDtweedie)

# load data set
data(auto)

# 5-fold cross validation using the lasso
cv0 <- cv.HDtweedie(x=auto$x,y=auto$y,p=1.5,nfolds=5,lambda.factor=.0005)

# make a CV plot
plot(cv0)

# define group index
group1 <- c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)

# 5-fold cross validation using the grouped lasso 
cv1 <- cv.HDtweedie(x=auto$x,y=auto$y,group=group1,p=1.5,nfolds=5,lambda.factor=.0005)

# make a CV plot
plot(cv1)
}
\keyword{models}
\keyword{regression}
