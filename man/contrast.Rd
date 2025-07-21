\name{contrast}
\alias{contrast}
\title{Counterfactual contrast from the parametric or sequential g-formula for continuous multiple time point interventions}

\description{Estimation of a contrast between counterfactual outcomes under different values of (continuous) interventions, or across different time points, using the parametric or sequential g-formula.}
\usage{
contrast(X, abar, nodes, contrastType = "difference", measure = mean,
         cond = NULL, cilevel = 0.95, ...)
}

\arguments{
  \item{X}{
An object of class \code{gformula} produced by \code{\link{gformula}}, with option \code{ret = TRUE}, or \code{sgf}. 
}
  \item{abar}{
Numeric vector or the string 'natural'. Specifies the intervention value(s) for the contrast. If two values are given, a contrast between these two intervention regimes is computed at the same outcome node. If a single value is given and \code{nodes} has two elements, a contrast between time points is computed under that intervention. If more than two entries are given, \code{contrastType} must be a custom function. See Details.
}
  \item{nodes}{
A character string vector specifying the variable(s) used in the contrast. If two values are given, a temporal contrast is computed (e.g., outcome change over time under the same intervention level). If more than two entries are given, \code{contrastType} must be a custom function. See Details.
}
  \item{contrastType}{
Type of contrast to compute between the counterfactual measures. Accepts one of 'difference', 'ratio', 'oddsratio', or a user-defined function taking \code{length(abar) * length(nodes)} numeric arguments and returning a numeric value. The default is 'difference'.
}
  \item{measure}{
Specifies the summary measure applied to the post-intervention counterfactual data. Defaults to \code{mean}.
}
  \item{cond}{
Optional filtering condition(s) applied to the post-intervention counterfactual data. Must be a quoted expression, e.g., \code{cond = quote(sex == 1)}, or a list of quoted expressions, e.g., \code{cond = list(quote(sex == 1), quote(sex == 0))}.
}
  \item{cilevel}{
Numeric value between 0 and 1 specifying the confidence level of the bootstrap confidence intervals. Defaults to 95\%. 
}
  \item{...}{
Additional arguments to be passed to \code{measure}.
}
}

\details{
Causal effects are defined as contrasts between the distributions of counterfactual variables under different interventions, across different time points or across different covariate strata. The counterfactual distributions to compare must be uniquely determined, by either specifying two values of \code{abar} at a single \code{nodes} or two \code{nodes} at a single intervention level \code{abar} or the natural course scenario with \code{abar = 'natural'} or two covariate strata \code{cond}. If the natural course scenario is selected and two \code{nodes} are specified, the natural intervention is compared across the two nodes. If one \code{nodes} is specified, the natural and observed scenarios are compared at a single node.

By default, the difference between the expectations of the two counterfactual outcome distributions is calculated. The difference can be exchanged for a ratio, odds ratio or custom contrast in the \code{contrastType} argument, and expectations can be exchanged for custom measures in the \code{measure} argument. Conditional measures can be specified through the \code{cond} argument. Custom contrasts, including those comparing more than two counterfactuals, can be defined by passing a function to \code{contrastType}.

Confidence intervals are based on the nonparametric bootstrap with \code{B} samples.
}

\value{
Returns a list of class \code{contrastResult}:
  \item{counterfactuals}{The estimated \code{measure}s of the counterfactual distributions.}
  \item{contrast}{The estimated contrast between the counterfactual measures.}
  \item{ciContrast}{The lower and upper bounds of the bootstrap confidence interval for the contrast.}
  \item{B}{The number of successful bootstrap samples. Will usually be equal to the input \code{B}.}
  \item{varContrast}{The estimated bootstrap variance of the contrast.}
}

\seealso{\code{\link{gformula}} and \code{sgf} for estimating expected counterfactual outcomes under multiple intervention values and \code{\link{custom.measure}} for measures other than expectations.}

\examples{
\dontrun{
data(EFV)
gf1 <- gformula(
  X = EFV, Anodes = c("efv.0", "efv.1", "efv.2", "efv.3", "efv.4"),
  Ynodes = c("VL.0", "VL.1", "VL.2", "VL.3", "VL.4"),
  Lnodes = c("adherence.1", "weight.1", "adherence.2", "weight.2",
             "adherence.3", "weight.3", "adherence.4", "weight.4"),
  abar = seq(1, 5), B = 10, ret = TRUE
)

# compare outcomes at last time point under (1,...,1) and (5,...,5) 
contrast(gf1, abar = c(1, 5), nodes = "VL.4")

# compare outcomes at different time points, for same intervention (2,...)
contrast(gf1, abar = 2, nodes = c("VL.3", "VL.2"))

# compare own measure (rel. risk reduction) instead of mean
# ... and conditional on subset
relativeRiskReduction <- function(k, l){(k - l) / k}

contrast(
  gf1, abar = c(1, 2), nodes = "VL.4",
  contrastType = relativeRiskReduction,
  cond = quote(sex == 1)
)

# Instead of the mean, any other measure can be taken,
# and - of course - applied also to counterfactual Lnodes
contrast(
  gf1, abar = 2, nodes = c("weight.3", "weight.2"),
  measure = median
)
}
}