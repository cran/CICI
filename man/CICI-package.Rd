\name{CICI-package}
\alias{CICI}
\docType{package}
\title{
Causal Inference with Continuous (Multiple Time Point) Interventions
}
\description{
This package facilitates the estimation of counterfactual outcomes for multiple values of continuous interventions at different time points, and allows plotting of causal dose-response curves. 

It implements the standard g-methods approach using the (semi-)parametric \emph{g}-formula, as described in the Schomaker et al. (2024) reference listed below. Weighted dose-response curves that address positivity violations, and are fitted via sequential \emph{g}-computation, are currently only available on GitHub and are not (yet) integrated in this package. 

The main function of the package is currently \code{\link{gformula}}.
}

\details{
\tabular{ll}{
Package: \tab CICI\cr
Type: \tab Package\cr
Version: \tab 0.9.4\cr
Date: \tab 2024-12-10\cr
License: \tab GPL-2\cr
Depends: \tab R (>= 4.0)\cr
Imports: \tab mgcv, glmnet, ggplot2, parallel, doParallel, foreach, doRNG, rngtools\cr
}

}
\author{
Michael Schomaker

Maintainer: Michael Schomaker <michael.schomaker@stat.uni-muenchen.de>
}

\references{
Schomaker M, McIlleron H, Denti P, Diaz I. (2024) \emph{Causal Inference for Continuous Multiple Time Point Interventions}, Statistics in Medicine, in press, see also \emph{https://arxiv.org/abs/2305.06645}.
}
