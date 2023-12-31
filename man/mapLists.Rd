\name{mapLists}
\alias{mapLists}
\title{Replaces contents of list A with elements of list B...}
\usage{mapLists(A, B, removeNAs=TRUE)}
\description{Replaces contents of list A with elements of list B}
\details{Combines two lists, \code{A} and \code{B}, such that \code{names(A)} are preserved, mapping to the
values of \code{B}, using \code{names(B)} as look up.  Ie. replaces values in \code{A} with values
in \code{B}, using \code{names(B)} as look up for values in \code{A}.
Once more?  See examples.
\emph{NB!} None-mapped entries are returned as NA, but can be removed using \code{\link{removeNAs}}.}
\value{List.}
\seealso{\code{\link{removeNAs}}}
\author{Stefan McKinnon Edwards \email{stefan.hoj-edwards@agrsci.dk}}
\arguments{\item{A}{List, elements are coerced to character for mapping to B.}
\item{B}{List.}
\item{removeNAs}{Boolean, whether to remove the \code{NA}s that occur because an element was not found in \code{B}.}}
\examples{A <- list('a1'='alpha','a2'='beta','a3'=c('gamma','delta'))
B <- list('alpha'='b1', 'gamma'=c('b2', 'b3'), 'delta'='b4')
mapLists(A, B)}
