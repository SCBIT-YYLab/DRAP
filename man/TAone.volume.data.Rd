\name{TAone.volume.data}
\docType{data}
\alias{TAone.volume.data}
\title{Example of volume data for TAone pattern}
\description{
An example of volume data for TAone pattern, which includes the required imformation for DRAP package.
}
\usage{
data(TAone.volume.data)
}
\details{
The dataset derived from Novartis Institutes for BioMedical Research PDX encyclopedia (NIBR PDXE), which includes tumor volume data for 6 tumor types, 277 tumors, and total 4771 animals responded to 61 treatments. The dataset used in DRAP includes information "Tumor" "ID" "Type" "Arms" "Times" "Volume".
}
\references{
Gao, H., et al. High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response. Nat Med 2015;21(11):1318-1325.
}
\examples{
data(TAone.volume.data)
TAone.volume.data[1:10,]
}
\keyword{datasets}