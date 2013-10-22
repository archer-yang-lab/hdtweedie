\name{print.HDtweedie}
\alias{print.HDtweedie}
\title{print a HDtweedie object}
\description{
Print the nonzero group counts at each lambda along the HDtweedie path.
}
\usage{
\method{print}{HDtweedie}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted \code{\link{HDtweedie}} object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
Print the information about the nonzero group counts at each lambda step in the \code{\link{HDtweedie}} object. The result is a two-column matrix with columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of the groups that have nonzero within-group coefficients, the \code{Lambda} column is the the corresponding lambda.
}
\value{
a two-column matrix, the first columns is the number of nonzero group counts and the second column is \code{Lambda}.}

\author{Wei Qian, Yi Yang and Hui Zou\cr
Maintainer: Wei Qian  <weiqian@stat.umn.edu>}
\references{
Qian, W., Yang, Y., Yang, Y. and Zou, H. (2013), ``Tweedie's Compound
    Poisson Model With Grouped Elastic Net,'' submitted to \emph{Journal of Computational and Graphical Statistics}.\cr
}
\examples{
# load HDtweedie library
library(HDtweedie)

# load auto data set
data(auto)

# fit the lasso
m0 <- HDtweedie(x=auto$x,y=auto$y,p=1.5)

# print out results
print(m0)

# define group index
group1 <- c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)

# fit the grouped lasso
m1 <- HDtweedie(x=auto$x,y=auto$y,group=group1,p=1.5)

# print out results
print(m1)
}
\keyword{models}
\keyword{regression}

