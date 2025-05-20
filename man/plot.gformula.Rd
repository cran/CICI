\name{plot.gformula}
\alias{plot.gformula}
\title{
Plot dose-response curves
}
\description{Function to plot dose-response curves based on results returned from \code{\link{gformula}}}
\usage{
\method{plot}{gformula}(x, msm.method = c("line","loess", "gam", "none"), 
                        CI = FALSE, time.points = NULL, 
                        cols = NULL, weight = NULL, xaxis=NULL,
                        variable = "psi", difference = FALSE, ...)
}
\arguments{
  \item{x}{
An object of \code{class} `gformula'.
}
  \item{msm.method}{
A string specifying the method to connect individual estimates into a curve (marginal structural model). One of \code{"line"},\code{"none"},\code{"gam"} and \code{"loess"}.}
  \item{CI}{
Logical. If TRUE, confidence bands are drawn; or confidence intervals for specific points if both \code{msm.method="none"} and appropriate.
}
  \item{time.points}{
A vector of time points for which the respective curves should be drawn. Default is all time points.
}
  \item{cols}{
A vector of strings specifying custom colours for each drawn curve.
}
  \item{weight}{
Weight vector of size "number of interventions times time points", that is used for the MSM if \code{msm.method="loess"} or \code{msm.method="gam"}.
}
\item{xaxis}{Either NULL or a string. If set to "time", then the x-axis is forced to represent time (unless this is impossible)}

  \item{variable}{
A string specifying the variable to be plotted under the natural course scenario (i.e., if \code{abar"natural"} and \code{cbar="natural"} in the respective \code{gformula} object). 
}
  \item{difference}{
Logical. If TRUE, differences of observed outcomes and outcomes under the natural intervention will be plotted (if \code{abar"natural"} and \code{cbar="natural"} in the respective \code{gformula} object.).
}
\item{...}{Further arguments to be passed on}
}

\details{
Time points and variable names should be specified according to the labeling of the results table returned by \code{\link{gformula}}.
}
\value{Draws an object of \code{class} `ggplot'.}

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


plot(est)
plot(est, time.points=c(1,5))
}
}

\keyword{dose-response curve}
