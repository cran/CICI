\name{model.formulas.update}
\alias{model.formulas.update}
\title{Update model formulas based on variable screening}
\description{Wrapper function to facilitate variable screening on all models generated through \code{\link{make.model.formulas}} and return updated formulas in the appropriate format for \code{\link{gformula}}. }
\usage{
model.formulas.update(formulas, X, screening = screen.glmnet.cramer,
                      with.s = FALSE, by= NA, ...)
}
\arguments{
  \item{formulas}{
A named list of length 4 containing model formulas for all Y-/L-/A- and Cnodes. These are likely formulas returned from \code{\link{make.model.formulas}}.
}
  \item{X}{
A data frame on which the model formulas are to be evaluated.}
  \item{screening}{
A screening function. Default is \code{screen.glmnet.cramer}, see Details below.}
  \item{with.s}{
Logical. If TRUE, a spline, i.e. s(), will be added to \emph{all} continuous variables.
}
\item{by}{
A character vector specifying the variables with which to multiply the smooth (if \code{with.s=TRUE}). 
}
  \item{\dots}{
optional arguments to be passed to the screening algorithm
}
}
\details{
The default screening algorithm uses LASSO for variable screening (and Cramer's V for the categorized version of all variables if LASSO fails). It is possible to provide user-specific screening algorithms. 
User-specific algorithms should take the data as first argument, \emph{one} model formula (i.e. one entry of the list in \code{model.formulas}) as second argument and return a vector of strings, containing the variable names that remain after screening. Another screening algorithm available in the package is \code{screen.cramersv}, which categorizes all variables, calculates their association with the outcome based on Cramer's \emph{V} and selects the 4 variables with strongest associations (can be changed with option \code{nscreen}).
The manual provides more information.

The fitted models of the updated models can be evaluated with \code{\link{fit.updated.formulas}}.
}
\value{
A list of length 4 containing the updated model formulas:
 \item{Lnames}{A vector of strings containing updated model formulas for all L nodes.}
 \item{Ynames}{A vector of strings containing updated model formulas for all Y nodes.}
 \item{Anames}{A vector of strings containing updated model formulas for all A nodes.}
 \item{Cnames}{A vector of strings containing updated model formulas for all C nodes.}
}

\seealso{
\code{\link{make.model.formulas}}, \code{\link{model.update}}, \code{\link{fit.updated.formulas}}
}

\examples{
\donttest{
data(EFV)

# first: generate generic model formulas
m <- make.model.formulas(X=EFV,
                         Lnodes  = c("adherence.1","weight.1",
                                     "adherence.2","weight.2",
                                     "adherence.3","weight.3",
                                     "adherence.4","weight.4"
                                    ),
                         Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                         Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                         evaluate=FALSE) 
                         
# second: update these model formulas based on variable screening with LASSO
glmnet.formulas <-  model.formulas.update(m$model.names, EFV)
glmnet.formulas 


# third: use these models for estimation
est <- gformula(X=EFV,
                Lnodes  = c("adherence.1","weight.1",
                            "adherence.2","weight.2",
                            "adherence.3","weight.3",
                            "adherence.4","weight.4"
                ),
                Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                Yform=glmnet.formulas$Ynames, Lform=glmnet.formulas$Lnames,
                abar=seq(0,2,1)
)
est
}
}

\keyword{model specification}