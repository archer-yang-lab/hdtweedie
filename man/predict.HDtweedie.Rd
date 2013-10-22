\name{predict.HDtweedie}
\alias{predict.HDtweedie}
\title{make predictions from a "HDtweedie" object.}

\description{
Similar to other predict methods, this functions predicts fitted values  from a \code{\link{HDtweedie}} object.
}

\usage{
\method{predict}{HDtweedie}(object, newx, s = NULL,
type=c("response","link"), ...)
}

\arguments{
		\item{object}{fitted \code{\link{HDtweedie}} model object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model. }
		\item{type}{
		type of prediction required:
		\itemize{
				\item Type \code{"response"} gives the mean response estimate.
				\item Type \code{"link"} gives the estimate for log mean response.}
		
		}
		\item{\dots}{Not used. Other arguments to predict. } }
\details{\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{predict} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \code{lambda} indices. }
\value{The object returned depends on type.}


\author{Wei Qian, Yi Yang and Hui Zou\cr
Maintainer: Wei Qian  <weiqian@stat.umn.edu>}
\references{
Qian, W., Yang, Y., Yang, Y. and Zou, H. (2013), ``Tweedie's Compound
    Poisson Model With Grouped Elastic Net,'' submitted to \emph{Journal of Computational and Graphical Statistics}.\cr
}
\seealso{\code{\link{coef}} method}
\examples{
# load HDtweedie library
library(HDtweedie)

# load auto data set
data(auto)

# fit the lasso
m0 <- HDtweedie(x=auto$x,y=auto$y,p=1.5)

# predicted mean response at x[10,]
print(predict(m0,type="response",newx=auto$x[10,]))

# define group index
group1 <- c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)

# fit the grouped lasso
m1 <- HDtweedie(x=auto$x,y=auto$y,group=group1,p=1.5)

# predicted the log mean response at x[1:5,]
print(predict(m1,type="link",newx=auto$x[1:5,]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
