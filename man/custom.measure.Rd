\name{custom.measure}
\alias{custom.measure}
\title{Custom estimands after applying \code{\link{gformula}}}
\description{The default estimate returned by \code{\link{gformula}} is the \strong{expected} outcome under the respective intervention strategies \code{abar}.
\code{custom.measure} takes an object of class \code{gformula} and enables estimation of other estimands based on the 
counterfactual datasets produced by \code{\link{gformula}} (if the option \code{ret=TRUE} had been chosen), for example estimands conditional on baseline variables, quantiles instead of expectations, and others.}
\usage{
custom.measure(X, fun = NULL, cond = NULL, verbose = TRUE, with.se = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
An object of class \code{gformula} produced by  \code{\link{gformula}} with option \code{ret=TRUE}.
}
  \item{fun}{
A function to be applied to the outcome(s) of the counterfactual data set.
}
  \item{cond}{
A string containing a condition to be applied to the counterfactual datasets.
}
  \item{verbose}{
Logical. TRUE if notes should be printed.
}
  \item{with.se}{Logical. TRUE if standard deviation should be calculated and returned.}
  \item{\dots}{
other parameters to be passed to \code{fun}
}
}
\details{In settings with censoring, it will often be needed to pass on the option \code{na.rm=T}, e.g. for the mean, median, quantiles, and others. 

Calculation of the bootstrap standard error (i.e., \code{with.se}=T) is typically not needed; but, for example, necessary for the calculations after multiple imputation and hence used by \code{\link{mi.boot}}.
}
\value{
An object of class \code{gformula}. See \code{\link{gformula}} for details.
}

\seealso{
see also \code{\link{gformula}}
}
\examples{
\donttest{
data(EFV)

est <- gformula(X=EFV,
                Lnodes  = c("adherence.1","weight.1",
                            "adherence.2","weight.2",
                            "adherence.3","weight.3",
                            "adherence.4","weight.4"
                ),
                Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                abar=seq(0,2,1), ret=TRUE
)

est
custom.measure(est, fun=prop,categ=1) # identical
custom.measure(est, fun=prop,categ=0)
custom.measure(est, fun=prop, categ=0, cond="sex==1")
# note: metabolic has been recoded internally (see output above)
custom.measure(est, fun=prop, categ=0, cond="metabolic==0") 
# does not make sense here, just for illustration (useful for metric outcomes)
custom.measure(est, fun=quantile, probs=0.1) 
}
}

% R documentation directory.
\keyword{estimands}
