\name{make.model.formulas}
\alias{make.model.formulas}
\title{
Compose appropriate model formulas
}
\description{
Function that generates generic model formulas for Y-/L-/A- and Cnodes, according to time ordering and to be used in \code{\link{gformula}} or \code{\link{model.formulas.update}}.
}
\usage{
make.model.formulas(X, Ynodes = NULL, Lnodes = NULL, Cnodes = NULL, Anodes = NULL,
                    survival = FALSE, evaluate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame, following the time-ordering of the nodes.  
}
  \item{Ynodes}{
A character string of column names in \code{X} of  the outcome variable(s).}
  \item{Lnodes}{
A character string of column names in \code{X} of  time-dependent (post first treatment) variable(s).
}
  \item{Cnodes}{
A character string of column names in \code{X} of  the censoring variable(s).
}
  \item{Anodes}{
A character string of column names in \code{X} of intervention variable(s).
}
  \item{survival}{
Logical. If TRUE, a survival setting is assumed and taken into account for model specification.
  }
  \item{evaluate}{
Logical. TRUE if model formulas should model formulas be evaluated on \code{X}.
}
}
\details{
This is a helper function to generate model formulas for Y-/L-/A- and Cnodes, according to the time ordering: i.e. to generate GLM/GAM model formulas for the respective nodes given all \emph{past} variables. In survival settings, past censoring and outcome nodes are omitted from the formulae. If censoring is present without a survival setting (e.g. Cnodes describe drop-outs and Y is a continuous outcome), then survival should be set as FALSE.
}
\value{Returns a named list:
\item{model.names}{A list of length 4 containing strings of the actual formulas}
\item{fitted.models}{A list of the fitted models (if \code{evaluate=TRUE})}
\item{fitted.model.summary}{A list of the summary of the fitted models (if \code{evaluate=TRUE})}
}

\seealso{The generated generic model formulas can be updated manually with \code{\link{model.update}} or in an automated manner with screening using \code{\link{model.formulas.update}}.}

\examples{
data(EFV)

m <- make.model.formulas(X=EFV,
                         Lnodes  = c("adherence.1","weight.1",
                                     "adherence.2","weight.2",
                                     "adherence.3","weight.3",
                                     "adherence.4","weight.4"
                                    ),
                         Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                         Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                         evaluate=FALSE) # set TRUE to see fitted models
                         
m$model.names # all models potentially relevant for gformula(), given full past
}

\keyword{model specification}