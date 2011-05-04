\name{AnnotationFuncs-package}
\title{Annotation translation functions}
\description{Annotation translation functions}
\details{\tabular{ll}{
Package:  \tab  AnnotationFuncs\cr
Type:  \tab  Package\cr
Version:  \tab  1.1.2 \cr
Date:  \tab  2011-05-04\cr
License:  \tab  Artistic (>=  2)\cr
LazyLoad:  \tab  yes\cr
}

Functions for handling translations between different identifieres using
the Biocore Data Team data-packages (e.g. \code{org.Bt.eg.db}).  
Primary function is \code{\link{translate}} for translating.  
Other functions include functions for selecting Refseqs or Gene Ontologies (GO).}
\alias{AnnotationFuncs}
\docType{package}
\author{Stefan McKinnon Edwards  \email{stefanm.edwards@agrsci.dk}}
\references{\url{http://www.iysik.com/}}
\keyword{package}
\seealso{\code{\link{translate}}}
\note{Requires user to deliver the annotation packages such as \code{org.Bt.egREFSEQ}.}
\examples{library(org.Bt.eg.db)
gene.symbols <- c('DRBP1','SERPINA1','FAKE','BLABLA')
# Find entrez identifiers of these genes.
eg <- translate(gene.symbols, org.Bt.egSYMBOL2EG)
# Note that not all symbols were translated.

# Go directly to Refseq identifiers.
refseq <- translate(gene.symbols, from=org.Bt.egSYMBOL2EG, to=org.Bt.egREFSEQ)
# Pick the proteins:
pickRefSeq(refseq, priorities=c('NP','XP'), reduce='all')}