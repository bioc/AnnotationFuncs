\name{getHomolog}
\alias{getHomolog}
\title{Performs quicker lookup in homologe data packages...}
\usage{getHomolog(values, mapping, genus, threshold=1, from, to, ...)}
\description{Performs quicker lookup in homologe data packages}
\details{Using the INPARANOID data packages such as \code{hom.Hs.inp.db} is very, very slow and can take up to 11 min (on the developers workstation).
This function introduces a new method that can do it in just 20 seconds (on the developers workstation).
In addition, it includes options for translating between different identifers both before and after the mapping.}
\value{List. Names of list corresponds to \code{values}, except those that could not be mapped nor translated.
Entries are character vectors.}
\references{\code{?hom.Hs.inp.db} - \url{http://inparanoid.sbc.su.se/}

Berglund, A.C., Sjolund, E., Ostlund, G., Sonnhammer, E.L.L. (2008) 
InParanoid 6: eukaryotic ortholog clusters with inparalogs
\emph{Nucleic Acids Res.} \bold{36}:D263--266

O'Brien, K.P., Maido, R., Sonnhammer, E.L.L (2005)
Inparanoid: A Comprehensive Database of Eukaryotic Orthologs
\emph{NAR} \bold{33}:D476--D480

Remm, M., Storm, C.E.V, Sonnhammer, E.L.L (2001)
Automatic clustering of orthologs and in-paralogs from pairwise species comparisons
\emph{J. Mol. Biol.} \bold{314}:1041--1052}
\seealso{\code{\link{translate}}, \code{\link{getTableName}}, \code{\link{mapLists}}}
\author{Stefan McKinnon Edwards \email{stefanm.edwards@agrsci.dk}}
\arguments{\item{values}{Vector, coerced to character vector, of values needed mapping by homology.}
\item{mapping}{Homology mapping object, such as \code{hom.Hs.inpBOSTA} or \code{revmap(hom.Hs.inpBOSTA)}.}
\item{genus}{Character vector. 5 character INPARANOID style genus name of the mapping object, e.g. 'BOSTA' for both \code{hom.Hs.inpBOSTA} and \code{revmap(hom.Hs.inpBOSTA)}.}
\item{thershold}{Numeric value between 0 and 1. Only clustered homologues with a parwise score above the threshold is included.
The native implementation has this set to 1.}
\item{pre.from}{Mapping object if \code{values} needs translation before mapping. 
E.g. \code{values} are entrez and \code{hom.Hs.inpBOSTA} requires ENSEMBLPROT, \code{hom.Hs.inpAPIME} requires Refseq (?).
Arguments \code{from} and \code{to} are just like in \code{\link{translate}}.}
\item{pre.to}{Second part of translation before mapping.}
\item{post.from}{Translate the result from homology mapping to a desired id; just like in \code{\link{translate}}.}
\item{post.to}{Second part of translation after mapping.}
\item{...}{Additional arguments sent to \code{\link{translate}}.}}
\examples{library(hom.Hs.inp.db)
library(org.Hs.eg.db)
library(org.Bt.eg.db)
getHomolog("ENSBTAP00000024572", revmap(hom.Hs.inpBOSTA), 'BOSTA') 
# And now, we will map from entrez genes 1, 2 and 3 to bovine Refseq
bovine.ensembl <- getHomolog(c(1,2,3), hom.Hs.inpBOSTA, 'BOSTA', from=org.Hs.egENSEMBLPROT, to=org.Bt.egENSEMBLPROT2EG)
refseqs <- translate(unlist(bovine.ensembl, use.names=F), org.Bt.egREFSEQ)
hs2bt.refseqs <- mapLists(bovine.ensembl, refseqs)
# Another way of doing it:
hs2bt.refseqs2 <- lapply(bovine.ensembl, translate, from=org.Bt.egREFSEQ, simplify=TRUE) # simplify=TRUE is very important here!}
