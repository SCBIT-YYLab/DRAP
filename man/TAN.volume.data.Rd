\name{TAN.volume.data}
\docType{data}
\alias{TAN.volume.data}
\title{Example of volume data for TAN pattern}
\description{
An example of volume data for TAN pattern, which includes the required imformation for DRAP package.
}
\usage{
data(TAN.volume.data)
}
\details{
Tumor volume data for TAN pattern must at least include the columns "Arms" "Tumor" "ID" "Times" and "Volume".
}
\examples{
data(TAN.volume.data)
TAN.volume.data[1:10,]
}
\keyword{datasets}