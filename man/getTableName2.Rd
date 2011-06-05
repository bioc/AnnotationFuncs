\name{.getTableName}
\alias{.getTableName}
\title{Gets the table name from the INPARANOID style genus names.}
\usage{.getTableName(genus)}
\description{Gets the table name from the INPARANOID style genus names.}
\details{The INPARANOID style genus name is a 5 letter acronym of the species name. 
Quote INPARANOID (\code{?hom.Hs.inpBOSTA}):

\emph{Names for these maps are done in the "INPARANOID style" which means that they are normally the 1st three letters of the genus followed by the 1st two letters of the species. For example: "Mus musculus" becomes "MUSMU", "Homo sapiens" becomes "HOMSA", "Monodelphis domestica" becomes "MONDO" etc. This means that for most of these organisms it will be possible to easily guess the abbreviations used. An exception may occur in the future if a new model organism has a very similar genus and species name to an existing one. }}
\value{Table name for genus.}
\references{\url{http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html}}
\author{Stefan McKinnon Edwards \email{stefanm.edwards@agrsci.dk}}
\arguments{\item{genus}{5 character INPARANOID genus name.}}
\examples{.getTableName('BOSTA')
.getTableName('mondo')}
