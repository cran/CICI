\name{EFVfull}
\alias{EFVfull}
\docType{data}
\title{
Pharmacoepidemiological HIV treatment data
}
\description{
A hypothetical, simulated dataset which is line with the data-generating process of Schomaker et al. (2024) and inspired by the data of Bienczak et al. (2017); see references below. Compared to the dataset \code{EFV}, it contains all variables of the DAG in Figure 3 of Schomaker et al. (2023), also those which are not needed for identification of the counterfactual quantity of interest; that is, the expected viral suppression (\code{VL}) under a specific intervention on efavirenz concentrations (\code{efv.0, efv.1, ...}).
}
\usage{data(EFVfull)}
\format{
  A data frame with 5000 observations on the following variables:
  \describe{
    \item{\code{sex}}{The patient's sex}
    \item{\code{metabolic}}{Metabolism status (slow, intermediate, extensive) related to the single nucleotide polymorphisms in the CYP2B6 gene, which is relevant
for metabolizing evafirenz and directly affects its concentration in the body.}
    \item{\code{log_age}}{log(age) at baseline}
    \item{\code{NRTI}}{Nucleoside reverse transcriptase inhibitor (NRTI) component of HIV treatment, i.e.
abacavir, stavudine or zidovudine.}
    \item{\code{weight.0}}{log(weight) at time 0 (baseline)}
    \item{\code{comorbidity.0}}{Presence of co-morbidities at time 0 (baseline)}
    \item{\code{dose.0}}{Dose of efavirenz administered at time 0 (basline)}
    \item{\code{efv.0}}{Efavirenz concentration at time 0 (baseline)}
    \item{\code{VL.0}}{Elevated viral load (viral failure) at time 0 (baseline)}
    \item{\code{adherence.1}}{Adherence at time 1 (if 0, then signs of non-adherence)}
    \item{\code{weight.1}}{log(weight) at time 1}
    \item{\code{comorbidity.1}}{Presence of co-morbidities at time 1}
    \item{\code{dose.1}}{Dose of efavirenz administered at time 1}
    \item{\code{efv.1}}{Efavirenz concentration at time 1}
    \item{\code{VL.1}}{Elevated viral load (viral failure) at time 1}
    \item{\code{adherence.2}}{Adherence at time 2 (if 0, then signs of non-adherence)}
    \item{\code{weight.2}}{log(weight) at time 2}
    \item{\code{comorbidity.2}}{Presence of co-morbidities at time 2}
    \item{\code{dose.2}}{Dose of efavirenz administered at time 2}
    \item{\code{efv.2}}{Efavirenz concentration at time 2}
    \item{\code{VL.2}}{Elevated viral load (viral failure) at time 2}
    \item{\code{adherence.3}}{Adherence at time 3 (if 0, then signs of non-adherence)}
    \item{\code{weight.3}}{log(weight) at time 3}
    \item{\code{comorbidity.3}}{Presence of co-morbidities at time 3}
    \item{\code{dose.3}}{Dose of efavirenz administered at time 3}
    \item{\code{efv.3}}{Efavirenz concentration at time 3}
    \item{\code{VL.3}}{Elevated viral load (viral failure) at time 3}
    \item{\code{adherence.4}}{Adherence at time 4 (if 0, then signs of non-adherence)}
    \item{\code{weight.4}}{log(weight) at time 4}
    \item{\code{comorbidity.4}}{Presence of co-morbidities at time 4}
    \item{\code{dose.4}}{Dose of efavirenz administered at time 4}
    \item{\code{efv.4}}{Efavirenz concentration at time 4}
    \item{\code{VL.4}}{Elevated viral load (viral failure) at time}
  }
}

\references{
Schomaker M, McIlleron H, Denti P, Diaz I. (2024) \emph{Causal Inference for Continuous Multiple Time Point Interventions}, Statistics in Medicine, 43:5380-5400, see also \emph{https://arxiv.org/abs/2305.06645}.

Bienczak et al. (2017) \emph{Determinants of virological outcome and adverse events in African children treated with paediatric nevirapine fixed-dose-combination tablets}, AIDS, \emph{31:905-915}
}
\examples{
data(EFVfull)
str(EFVfull)
}
\keyword{Datasets}
\keyword{efavirenz}