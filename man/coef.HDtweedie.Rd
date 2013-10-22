\name{coef.HDtweedie}
\alias{coef.HDtweedie}
\title{get coefficients or make coefficient predictions from an "HDtweedie" object.}
\description{
Computes the coefficients at the requested values for \code{lambda} from a fitted \code{\link{HDtweedie}} object.
}
\usage{
\method{coef}{HDtweedie}(object, s = NULL, ...)
}
\arguments{
		\item{object}{fitted \code{\link{HDtweedie}} model object.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which predictions are required. Default is the entire sequence used to create the model.}
		\item{\dots}{not used. Other arguments to predict. } 
}

\details{
\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{coef} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of coefficients from both left and right \code{lambda} indices. 
}
\value{The coefficients at the requested values for \code{lambda}.}


\author{Wei Qian, Yi Yang and Hui Zou\cr
Maintainer: Wei Qian  <weiqian@stat.umn.edu>}
\references{
Qian, W., Yang, Y., Yang, Y. and Zou, H. (2013), ``Tweedie's Compound
    Poisson Model With Grouped Elastic Net,'' submitted to \emph{Journal of Computational and Graphical Statistics}.\cr
}
\seealso{\code{\link{predict.HDtweedie}} method}
\examples{
# load HDtweedie library
library(HDtweedie)

# load data set
data(auto)

# fit the lasso
m0 <- HDtweedie(x=auto$x,y=auto$y,p=1.5)

# the coefficients at lambda = 0.01
coef(m0,s=0.01)

# define group index
group1 <- c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)

# fit grouped lasso
m1 <- HDtweedie(x=auto$x,y=auto$y,group=group1,p=1.5)

# the coefficients at lambda = 0.01 and 0.04
coef(m1,s=c(0.01,0.04))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
