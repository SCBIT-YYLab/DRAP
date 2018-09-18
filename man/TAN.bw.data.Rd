\name{TAN.bw.data}
\docType{data}
\alias{TAN.bw.data}
\title{Example of body weight data for TAN pattern}
\description{
An example of body weight data for TAN pattern, which includes the required imformation for DRAP package.
}
\usage{
data(TAN.bw.data)
}
\details{
Body weigth data for TAN pattern must at least include the columns "Arms" "Tumor" "ID" "Times" and "BodyWeight".
}
\examples{
data(TAN.bw.data)
TAN.bw.data[1:10,]
}
\keyword{datasets}