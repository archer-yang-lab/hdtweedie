\name{predict.cv.HDtweedie}
\alias{predict.cv.HDtweedie}
\title{make predictions from a "cv.HDtweedie" object.}
\description{
This function makes predictions from a cross-validated \code{HDtweedie} model,
using the stored \code{"cv.HDtweedie"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{predict}{cv.HDtweedie}(object, newx, s=c("lambda.1se","lambda.min"),...)
}
\arguments{
		\item{object}{fitted \code{\link{cv.HDtweedie}} object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix. See documentation for \code{predict.HDtweedie}.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the value \code{s="lambda.1se"} stored
		on the CV object. Alternatively \code{s="lambda.min"} can be
		used. If \code{s} is numeric, it is taken as the value(s) of
		\code{lambda} to be used.}
		\item{\dots}{not used. Other arguments to predict. } }
\details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The returned object depends on the \dots argument which is passed on
to the \code{\link{predict}} method for \code{\link{HDtweedie}} objects.}

\author{Wei Qian, Yi Yang and Hui Zou\cr
Maintainer: Wei Qian  <weiqian@stat.umn.edu>}
\references{
Qian, W., Yang, Y., Yang, Y. and Zou, H. (2013), ``Tweedie's Compound
    Poisson Model With Grouped Elastic Net,'' submitted to \emph{Journal of Computational and Graphical Statistics}.\cr
}

\seealso{\code{\link{cv.HDtweedie}}, and \code{\link{coef.cv.HDtweedie}} methods.}
\examples{
# load HDtweedie library
library(HDtweedie)

# load data set
data(auto)

# 5-fold cross validation using the lasso
cv0 <- cv.HDtweedie(x=auto$x,y=auto$y,p=1.5,nfolds=5) 

# predicted mean response at lambda = lambda.1se, newx = x[1,]
pre = predict(cv0, newx = auto$x[1,], type = "response")

# define group index
group1 <- c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)

# 5-fold cross validation using the grouped lasso 
cv1 <- cv.HDtweedie(x=auto$x,y=auto$y,group=group1,p=1.5,nfolds=5)

# predicted the log mean response at lambda = lambda.min, x[1:5,]
pre = predict(cv1, newx = auto$x[1:5,], s = cv1$lambda.min, type = "link")
}
\keyword{models}
\keyword{regression}
