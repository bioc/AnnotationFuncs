\name{removeNAs}
\alias{removeNAs}
\title{Removes entries equal NA from list or vector...}
\usage{removeNAs(l)}
\description{Removes entries equal \code{NA} from list or vector}
\details{Removes entries equal \code{NA}, but not mixed entries containing, amongst others, \code{NA}.  
Good for use after \code{\link{mapLists}} that might return entries equal \code{NA}.}
\author{Stefan McKinnon Edwards \email{stefanm.edwards@agrsci.dk}}
\arguments{\item{l}{Vector or list.}}
