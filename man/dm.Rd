\name{dm}
\alias{dm}
\title{ Fitting dynamic models }
\description{ dm is used to fit dynamic models.}
\usage{ dm( formula , data = parent.frame() , opts = NULL) }
\arguments{
	\item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.}
	\item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called}
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

garch11

}
\keyword{timeseries}
