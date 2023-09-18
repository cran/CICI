\name{mi.boot}
\alias{mi.boot}
\title{Obtaining estimates from multiply imputed data}
\description{Combines \code{\link{gformula}} estimates obtained from multiple imputed data sets according to the \emph{MI Boot} and \emph{MI Boot pooled} methods decribed in Schomaker and Heumann (2018, see reference section below) }
\usage{
mi.boot(x, fun, cond=NULL, pooled=FALSE, ...)
}
\arguments{
\item{x}{A list of objects of \code{class} `gformula'}
\item{fun}{A function to be applied to the outcome(s) of the counterfactual data set. For expected outcome, use \code{\link{mean}} and possibly pass on option \code{na.rm=TRUE}.}
\item{cond}{A string containing a condition to be applied to the counterfactual datasets.}
\item{pooled}{Logical. If TRUE, confidence interval estimation is based on the MI Boot pooled from Schomaker and Heumann (2018), otherwise on MI Boot.}
\item{...}{additional arguments to be passed on to \code{fun}}
}

\value{
An object of class \code{gformula}. See \code{\link{gformula}} for details.
}

\references{
Schomaker, M., Heumann, C. (2018) \emph{Bootstrap inference when using multiple imputation},
Statistics in Medicine, 37:2252-2266
}
\author{
Michael Schomaker
}

\examples{
\donttest{
data(EFV)

# suppose the following subsets were actually multiply imputed data (M=2)
EFV_1 <- EFV[1:2500,]
EFV_2 <- EFV[2501:5000,]

# first: conduct analysis on each imputed data set. Set ret=T.
m1 <- gformula(X=EFV_1,
               Lnodes  = c("adherence.1","weight.1",
                           "adherence.2","weight.2",
                           "adherence.3","weight.3",
                           "adherence.4","weight.4"
               ),
               Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
               Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
               abar=seq(0,5,1), verbose=FALSE, ret=TRUE
        )

m2 <- gformula(X=EFV_2,
               Lnodes  = c("adherence.1","weight.1",
                           "adherence.2","weight.2",
                           "adherence.3","weight.3",
                           "adherence.4","weight.4"
               ),
               Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
               Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
               abar=seq(0,5,1), verbose=FALSE, ret=TRUE
)

# second combine results
m_imp <- mi.boot(list(m1,m2), mean) # uses MI rules & returns 'gformula' object
plot(m_imp)

# custom estimand: evaluate probability of suppression (Y=0), among females
m_imp2 <- mi.boot(list(m1,m2), prop, categ=0, cond="sex==1")
plot(m_imp2)
}
}

% R documentation directory.
\keyword{multiple imputation}