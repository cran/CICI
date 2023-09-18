\name{model.update}
\alias{model.update}
\title{Update GAM models}
\description{A wrapper to simplify the update of GAM models}
\usage{
model.update(gam.object, form)
}
\arguments{
  \item{gam.object}{
A \code{gam} object produced with package \pkg{mgcv}. 
}
  \item{form}{
A new model formula in the form \code{.~formula}}
}
\details{The \code{gam} object needs to be fitted with the option control=list(keepData=T), otherwise the function can not access the data that is needed to update the model fit. Note that both \code{\link{fit.updated.formulas}} and \code{\link{make.model.formulas}} with option \code{evaluate=T} produce results that are based on this option.}
\value{An object of \code{class} `gam', `glm' and `lm'.
}

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
                         evaluate=TRUE) # set TRUE for model.update()

# update first confounder model of weight manually
model.update(m$fitted.models$fitted.L$m_weight.1, .~s(weight.0, by=sex))

# manual update of model formula
m$model.names$Lnames[2] <- "weight.1 ~ s(weight.0, by=sex)" 

}

\keyword{model specification}