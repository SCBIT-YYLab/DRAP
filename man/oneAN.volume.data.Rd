\name{oneAN.volume.data}
\docType{data}
\alias{oneAN.volume.data}
\title{Example of volume data for oneAN pattern}
\description{
An example of volume data for oneAN pattern, which includes the required imformation for DRAP package.
}
\usage{
data(oneAN.volume.data)
}
\details{
Tumor volume data for oneAN pattern must at least include the columns "Arms" "ID" "Times" and "Volume".
}
\examples{
data(oneAN.volume.data)
oneAN.volume.data[1:10,]
}
\keyword{datasets}