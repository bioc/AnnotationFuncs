\name{translate}
\alias{translate}
\title{Translate between different identifiers...}
\usage{translate(values, from, to, reduce=c("all", "first", "last"),
    return.list=TRUE, remove.missing=TRUE, simplify=FALSE, ...)}
\description{Translate between different identifiers}
\details{Function for translating from one annotation to another, eg. from RefSeq to  
Ensemble.  This function takes a vector of annotation values and translates     
first to the primary annotation in the Biocore Data Team package (ie. entrez gene identifier for org.Bt.eg.db)
and then to the desired product, while removing non-translated annotations   
and optionally reducing the result so there is only a one-to-one relation. 

If you want to do some further mapping on the result, you will have to use
either \code{unlist} og \code{lapply}, where the first returns all the end-products
of the first mapping, returning a new list, and the latter produces a list-within-list.

If \code{from} returns GO identifiers (e.g. \code{from = org.Bt.egGO}), then the 
returned resultset is more complex and consists of several layers of lists instead of 
the usual list of character vectors. If \code{to} has also been specified, the GO IDs 
must be extracted (internally) and you have the option of filtering for evidence and category at this point.
See \code{\link{pickGO}}.}
\note{Requires user to deliver the annotation packages such as org.Bt.egREFSEQ.}
\value{List; names of elements are \code{values} and the elements are the translated elements,
or \code{NULL} if not translatable with \code{remove.missing = TRUE}.}
\author{Stefan McKinnon Edwards \email{stefan.hoj-edwards@agrsci.dk}}
\seealso{\code{\link{pickRefSeq}}, \code{\link{pickGO}}}
\arguments{\item{values}{Vector of annotations that needs translation. Coerced to character vector.}
\item{from}{Type of annotation \code{values} are given in. NB! take care in the       
orientation of the package, ie. if you have RefSeq annotations, use  
\code{org.Bt.egREFSEQ2EG} or (in some cases) \code{revmap(org.Bt.egREFSEQ)}.}
\item{to}{Desired goal, eg. \code{org.Bt.egENSEMBLPROT}. If \code{NULL} (default), goal 
if the packages primary annotation (eg. entrez gene for org.Bt.eg.db).
Throws a warning if the organisms in \code{from} and \code{to} are not the same.}
\item{reduce}{Reducing method, either return all annotations (one-to-many relation)
or the first or last found annotation. The reducing step is applied  
after translating to the goal:                                       
\code{all}: returns all annotations                                       
\code{first} or \code{last}: choose first or last of arbitrarily ordered list.}
\item{return.list}{Logical, when \code{TRUE}, returns the translation as a list where names}
\item{remove.missing}{Logical, whether to remove non-translated values, defaults \code{TRUE}.}
\item{simplify}{Logical, unlists the result. Defaults to FALSE. Usefull when using \code{translate} in 
a \code{lapply} or \code{sapply}.}
\item{...}{Additional arguments sent to \code{\link{pickGO}} if \code{from} returns GO set.}}
\examples{library(org.Bt.eg.db)
genes <- c(280705, 280706, 100327208)
translate(genes, org.Bt.egSYMBOL)

symbols <- c("SERPINA1","KERA","CD5")
refseq <- translate(symbols, from=org.Bt.egSYMBOL2EG, to=org.Bt.egREFSEQ)
# Pick the proteins:
pickRefSeq(refseq, priorities=c('NP','XP'), reduce='all')

# If you wanted do do some further mapping on the result from 
# translate, simply use lapply.

library(GO.db)
GO <- translate(genes, org.Bt.egGO)
# Get all biological processes:
pickGO(GO, category='BP')
# Get all ontologies with experimental evidence:
pickGO(GO, evidence=c('IMP','IGI','IPI','ISS','IDA','IEP','IEA'))}
