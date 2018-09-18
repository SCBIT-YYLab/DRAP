\name{TAone.bw.data}
\docType{data}
\alias{TAone.bw.data}
\title{Example of body weight data for TAone pattern}
\description{
An example of body weight data for TAone pattern, which includes the required imformation for DRAP package.
}
\usage{
data(TAone.bw.data)
}
\details{
The dataset derived from Novartis Institutes for BioMedical Research PDX encyclopedia (NIBR PDXE), which includes tumor volume data for 6 tumor types, 277 tumors, and total 4771 animals responded to 61 treatments. The dataset used in DRAP includes information "Tumor" "ID" "Type" "Arms" "Times" "BodyWeight".
}
\references{
Gao, H., et al. High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response. Nat Med 2015;21(11):1318-1325.
}
\examples{
data(TAone.bw.data)
TAone.bw.data[1:10,]
}
\keyword{datasets}