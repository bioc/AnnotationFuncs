\name{.dbEscapeString}
\alias{.dbEscapeString}
\title{Private Escape string...}
\usage{.dbEscapeString(str, raise.error=TRUE)}
\description{Private Escape string}
\details{Does not escape strings, but raises an error if any character expect normal letters and underscores are found in the string.}
\value{Invisible logical}
\arguments{\item{str}{String to test}
\item{raise.error}{Logical, whether to raise an error or not.}}
