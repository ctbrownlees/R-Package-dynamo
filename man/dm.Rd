\name{dm}
\alias{dm}
\alias{print.dm}
\alias{plot.dm}
\title{Dynamic Modeling}
\description{dm}
\usage{ dm( formula , data = parent.frame(), param = NULL, opts = NULL) }
\arguments{
	\item{formula}{formula}
	\item{data}{data}
	\item{param}{param}
	\item{opts}{opts}
}
\details{The dm procedure can be used to fit dynamic models.}
\value{
	The 'dm' procedure returns an object of class 'dm'

	The functions 'summary', 'plot' and 'print' are used to obtain and print a summary of the results.

	An object of class "dm" is a list containing at least the following components:
}
\references{}
\author{Christian Brownlees}
\examples{

data( sp500 )

ret <- diff( log( sp500$price ) )*100

garch11 <- dm( ret ~ garch )

}
\keyword{timeseries}
