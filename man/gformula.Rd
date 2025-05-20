\name{gformula}
\alias{gformula}
\title{
Parametric g-formula for continuous multiple time point interventions
}
\description{
Estimation of counterfactual outcomes for multiple values of continuous interventions at different time points using the g-formula. 

}
\usage{
gformula(X, Anodes, Ynodes, Lnodes = NULL, Cnodes = NULL,
         abar = NULL, cbar = "uncensored",
         survivalY = FALSE, 
         Yform = "GLM", Lform = "GLM", Aform = "GLM", Cform = "GLM",
         calc.support = FALSE, B = 0, ret = FALSE, ncores = 1, 
         verbose = TRUE, seed = NULL, prog = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame, following the time-ordering of the nodes. Categorical variables with k categories should be a factor, with levels 0,...,k-1. Binary variables should be coded 0/1. 
}
  \item{Anodes}{
A character string of column names in \code{X} of the intervention variable(s).
}
  \item{Ynodes}{
A character string of column names in \code{X} of  the outcome variable(s).
}
  \item{Lnodes}{
A character string of column names in \code{X} of  the time-dependent (post first treatment) variable(s).
}
  \item{Cnodes}{
A character string of column names in \code{X} of  the censoring variable(s).
}
  \item{abar}{
Numeric vector or matrix of intervention values, or the string "natural". See Details.
}
  \item{cbar}{
Typically either the string "uncensored" or "natural", but a numeric vector or matrix of censoring values is not forbidden. See Details.
}
  \item{survivalY}{
Logical. If TRUE, then Y nodes are indicators of an event, and if Y at some time point is 1, then all following should be 1.
}
  \item{Yform}{
A string of either "GLM", "GAM" or of length 'number of Ynodes' with model formulas. See Details. 
}
  \item{Lform}{
A string of either "GLM", "GAM" or of length 'number of Lnodes' with model formulas. See Details. 
}
  \item{Aform}{
A string of either "GLM", "GAM" or of length 'number of Anodes' with model formulas. See Details. 
}
  \item{Cform}{
A string of either "GLM", "GAM" or of length 'number of Cnodes' with model formulas. See Details. 
}
  \item{calc.support}{
Logical. If \code{TRUE}, both crude and conditional support is estimated.
}
  \item{B}{
An integer specifying the number of bootstrap samples to be used, if any. 
}
  \item{ret}{
Logical. If \code{TRUE}, the simulated post-intervention data is returned.
}
  \item{ncores}{
An integer for the number of threads/cores to be used. If >1, parallelization will be utilized.
}
  \item{verbose}{
Logical. If \code{TRUE}, notes and warnings are printed.
}
  \item{seed}{
An integer specifying the seed to be used to create reproducable results for parallel computing (i.e. when ncores>1).
}
  \item{prog}{
A character specifying a path where progress should be saved (typically, when \code{ncores}>1)
}
  \item{...}{
Further arguments to be passed on.
}
}
\details{
By default, expected counterfactual outcomes (specified under \code{Ynodes}) under the intervention \code{abar} are calculated. Other estimands can be specified via \code{\link{custom.measure}}.

If \code{abar} is a vector, then each vector component is used as the intervention value at each time point; that is, interventions which are constant over time are defined. If \code{abar} is a matrix (of size 'number interventions' x 'time points'), then each row of the length of \code{Anodes} refers to a particular time-varying intervention strategy. The natural intervention can be picked by setting \code{abar='natural'}.

The fitted outcome and confounder models are based on generalized additive models (GAMs) as implemented in the \code{mgcv} package. Model families are picked automatically and reported in the output if \code{verbose=TRUE} (see manual for modifications, though they hardly ever make sense). The model formulas are standard GLMs or GAMs (with penalized splines for continuous covariates), conditional on the past, unless specific formulae are given. It is recommended to use customized formulae to reduce the risk of model mis-specification and to ensure that the models make sense (e.g., not too many splines are used when this is computationally not meaningful). This can be best facilitated by using objects generated through \code{\link{make.model.formulas}}, followed by \code{\link{model.formulas.update}} and/or \code{\link{model.update}} (see examples for those functions). 

For survival settings, it is required that i) \code{survivalY=TRUE} and ii) after a Cnode/Ynode is 1, every variable thereafter is set to \code{NA}. See manual for an example. By default, the package intervenes on Cnodes, i.e. calculates counterfactual outcomes under no censoring. 

If \code{calc.support=TRUE}, conditional and crude support measures (i.e., diagnostics) are calculated as described in Section 3.4 of Schomaker et al. (2023). Another useful diagnostic for multiple time points is the natural course scenario, which can be evaluated under \code{abar='natural'} and \code{cbar='natural'}.

To parallelize computations automatically, it is sufficient to set \code{ncores>1}, as appropriate. To make estimates under parallelization reproducible, use the \code{seed} argument. To watch the progress of parallelized computations, set a path in the \code{prog} argument: then, a text file reports on the progress, which is particularly useful if lengthy bootstrapping computations are required.
}
\value{
Returns an object of of \code{class} `gformula':
  \item{results}{data.frame of results. That is, the estimated counterfactual outcomes depending on the chosen intervention strategies, and time points.}
  \item{diagnostics}{list of diagnostics and weights based on the estimated support (if \code{calc.support=TRUE})}
  \item{simulated.data}{list of counterfactual data sets related to the interventions defined through option \code{abar} (and \code{cbar}). Will be \code{NULL} is \code{ret=FALSE}.}
  \item{observed.data}{list of observed data (and bootstrapped observed data). Will be \code{NULL} is \code{ret=FALSE}.}
  \item{setup}{list of chosen setup parameters}
}

\author{
Michael Schomaker
}


\seealso{
\code{\link{plot.gformula}} for plotting results as (causal) dose response curves, \code{\link{custom.measure}} for evaluating custom estimands and \code{\link{mi.boot}} for using \code{\link{gformula}} on multiply imputed data.}

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
                        abar=seq(0,10,1)
)
est
}

}
\keyword{g-formula}
