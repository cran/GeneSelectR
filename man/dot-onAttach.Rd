% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{.onAttach}
\alias{.onAttach}
\title{Package Attachment Function}
\usage{
.onAttach(libname, pkgname)
}
\arguments{
\item{libname}{The name of the library.}

\item{pkgname}{The name of the package.}
}
\description{
This function is called when the package is attached. It checks for the availability of essential Bioconductor packages and alerts if any are missing.
}
\details{
The function checks for the presence of specific Bioconductor packages that are essential for the package's functionality. If any required packages are missing, it displays a startup message advising the user to install the missing packages using \code{BiocManager::install()}.

Instead of stopping the package loading process, it alerts the user about any missing dependencies, recommending their installation for full functionality.
}
\keyword{internal}
