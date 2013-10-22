\name{plot.HDtweedie}
\alias{plot.HDtweedie}
\title{Plot solution paths from a "HDtweedie" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{HDtweedie}} object.
}
\usage{
\method{plot}{HDtweedie}(x, group = FALSE, log.l = TRUE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{HDtweedie}} model}
		\item{group}{what is on the Y-axis. Plot the norm of each group
		 if \code{TRUE}. Plot each coefficient if \code{FALSE}.}
		\item{log.l}{what is on the X-axis. Plot against the log-lambda
		sequence if \code{TRUE}. Plot against the lambda sequence if \code{FALSE}.}
		\item{\dots}{other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced.
}
\author{Wei Qian, Yi Yang and Hui Zou\cr
Maintainer: Wei Qian  <weiqian@stat.umn.edu>}
\references{
Qian, W., Yang, Y., Yang, Y. and Zou, H. (2013), ``Tweedie's Compound
    Poisson Model With Grouped Elastic Net,'' submitted to \emph{Journal of Computational and Graphical Statistics}.\cr
}
\examples{
# load HDtweedie library
library(HDtweedie)

# load data set
data(auto)

# fit the lasso
m0 <- HDtweedie(x=auto$x,y=auto$y,p=1.5)

# make plot
plot(m0) # plots the coefficients against the log-lambda sequence

# define group index
group1 <- c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)

# fit group lasso
m1 <- HDtweedie(x=auto$x,y=auto$y,group=group1,p=1.5)

# make plots
par(mfrow=c(1,3))
plot(m1) # plots the coefficients against the log-lambda sequence 
plot(m1,group=TRUE) # plots group norm against the log-lambda sequence 
plot(m1,log.l=FALSE) # plots against the lambda sequence
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
