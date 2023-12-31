\name{fit.updated.formulas}
\alias{fit.updated.formulas}
\title{Fit models after screening}
\description{Fits the models that have been generated with screening using \code{\link{model.formulas.update}}.}
\usage{
fit.updated.formulas(formulas, X)
}
\arguments{
  \item{formulas}{An object returned by \code{\link{model.formulas.update}}}
  \item{X}{A data frame on which the model formulas should be evaluated}
}
\details{
Fits generalized (additive) linear models based on the screened model formula list generated by \code{\link{model.formulas.update}}.
}
\value{
Returns a list of length 2:
\item{fitted.models}{A list of length 4, containing the fitted Y-/L-/C- and A-models.}
\item{all.summaries}{A list of length 4, containing the summary of the fitted Y-/L-/C- and A-models.}
}

\seealso{\code{\link{model.formulas.update}}}

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

# then: fit and inspect the updated models
fitted.models <- fit.updated.formulas(glmnet.formulas, EFV)
fitted.models$all.summaries
fitted.models$all.summaries$Ynames[1] # first outcome model
}
}

\keyword{model specification}
