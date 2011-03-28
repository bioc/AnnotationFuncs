\name{pickGO}
\alias{pickGO}
\title{Cleans up result from org...}
\usage{pickGO(l, evidence=NA, category=NA)}
\description{Cleans up result from org.Xx.egGO and returns specific GO identifiers}
\details{Cleans up result from org.Xx.egGO and returns GO identifier for  either
biological process (BP), cellular component (CC), or molecular function (MF).  
Can be used on list of GOs from \code{\link{translate}}, or a single list of GOs from an annotation package.  
May reduce list, if the (sub)list does not contain the chosen class!}
\value{List with only the picked elements.}
\author{Stefan McKinnon Edwards \email{stefanm.edwards@agrsci.dk}}
\seealso{\code{\link{pickRefSeq}}, \code{\link{getEvidenceCodes}}, \code{\link{translate}}}
\arguments{\item{l}{Character vector, or list of, og GO identifiers.}
\item{evidence}{Character vector, filters on which kind of evidence to return; for a larger list see \code{\link{getEvidenceCodes}}. \\*
Evidence codes may be: \code{c('IMP','IGI','IPI','ISS','IDA','IEP','IEA','TAS','NAS','ND','IC')}. \\*
Leave as \code{NA} to ignore filtering on this part.}
\item{category}{Character vector, filters on which ontology to return: biological process (BP), cellular component (CC), or molecular function (MF). \\*
Leave as \code{NA} to ignore filtering on this part.}}
\examples{library(org.Bt.eg.db)
genes <- c(280705, 280706, 100327208)
GO <- translate(genes, org.Bt.egGO)
# Get all biological processes:
pickGO(GO, category='BP')
# Get all ontologies with experimental evidence:
pickGO(GO, evidence=c('IMP','IGI','IPI','ISS','IDA','IEP','IEA'))}
