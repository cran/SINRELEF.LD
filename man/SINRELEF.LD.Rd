\name{SINRELEF.LD}
\alias{SINRELEF.LD}

\title{
Score Information, reliability, and relative efficiency in measures that contain locally-dependent items.
}
\description{
\code{SINRELEF.LD} implements an approach aimed at assessing the accuracy and effectiveness of raw scores obtained in scales that contain locally dependent items. The program uses as input the calibration (structural) item estimates obtained from fitting extended unidimensional factor-analytic solutions in which the existing local dependencies are included. Measures of reliability (Omega) and information are proposed at three levels: (a) total score, (b) bivariate-doublet, and (c) item-by-item deletion, and  are compared to those that would be obtained if all the items had been locally independent. All the implemented procedures can be obtained from: (a) linear factor-analytic solutions in which the item scores are treated as approximately continuous, and (b) non-linear solutions in which the item scores are treated as ordered-categorical.
}
\usage{

SINRELEF.LD(L, PSI, THRES, ncat, model = 'linear', doublet_list, cor_doublet, N,
  CI = 90, display = TRUE)
}
\arguments{
  \item{L}{
A vector containing the item loading estimates in the unidimensional factor-analytic (FA) solution.
}
  \item{PSI}{
A vector containing the item residual standard deviations.
}
  \item{THRES}{
A vector containing the item thresholds.
}
  \item{ncat}{
Number of item response categories (for the graded model only).
}
  \item{model}{
The factor-analytic model used for calibrating the item set, being 'linear' for the standard FA model or 'graded' for the non-linear ordered-categorical FA model.
}
  \item{doublet_list}{
A vector containing the pairs of items with correlated specificities.
}
  \item{cor_doublet}{
A vector containing the residual correlations corresponding to each pair, specified in doublet_list.
}
  \item{N}{
The number of observations in the sample.
}
  \item{CI}{
Choice of the confidence level for the intervals (90 or 95, 90 by default).
}
  \item{display}{
Determines if the output will be displayed in the console, TRUE by default. If it is TRUE, the output is returned silently and if it is FALSE, the output is returned in the console.
}
}
\details{

\code{SINRELEF.LD} implements a series of procedures developed by the authors of the package, which have two main purposes: the first is to assess the real amount of reliability and information provided by the item and scale scores when the scale contains locally dependent items. The second is to compare the results of the assessment to the results that would have been obtained if all the items had been locally independent.

The procedures are intended at three score levels. First, at the total score, \code{SINRELEF.LD} provides the Omega reliability estimate and the amount of score information obtained taking in to account the local dependencies as well as the estimates obtained if local independence had been assumed. Furthermore, a relative efficiency index that assesses the loss of efficiency of the scores due to the local dependencies is provided. Second, at the bivariate level, for each specified doublet, \code{SINRELEF.LD} provides the corresponding relative efficiency estimate. Finally, at the single item level, \code{SINRELEF.LD} provides estimates of the change in reliability and relative information if the item was omitted from the scale.

The indices implemented in the function assume that the items have been calibrated using an extended factor-analytic solution that models the existing local dependencies. So, the input of the function requieres the calibration estimates obtained by fitting a solution of this type (e.g. loadings, residual variances, and residual correlations) to be provided. This type of extended solutions can be obtained from widely available programs such as \code{Mplus}, \code{FACTOR} or \code{lavaan} package. Furthermore, the proposed approach works with solutions obtained from two types of factor-analytic models. First, is the linear model, in which the item scores are treated as approximately continuous. Second, is the non-linear model, based on the underlying variables approach in which the item scores are treated as orderered-categorical.

Finally, confidence intervals based on the Bollen-Stine empirical approach are proposed for all the indices provided in the package.
}
\value{
\item{omld}{Omega reliability estimate taking in to account the item local dependences.}
\item{omli}{Omega reliability estimate if all the items were locally independent.}
\item{relef}{Score relative efficiency.}
\item{relef_doublet}{Doublet score relative efficiency.}
\item{omega_del}{Reliability estimate if the item was omitted.}
\item{r_info_del}{Relative information if the item was omitted.}

}
\references{

Lorenzo-Seva, U., & Ferrando, P. J. (2006). FACTOR: A computer program to fit the exploratory factor analysis model. Behavior research methods, 38(1), 88-91. \doi{10.3758/BF03192753}

Muthén, L. K., & Muthén, B. O. (1998-2018). Mplus User's Guide. Sixth Edition. Los Angeles, CA: Muthén & Muthén.

Rosseel Y (2012). “lavaan: An R Package for Structural Equation Modeling.” Journal of Statistical Software, 48(2), 1–36.\doi{10.18637/jss.v048.i02}.
}
\author{
David Navarro-Gonzalez

Fabia Morales-Vives

Pere Joan Ferrando
}
\examples{
## Example data set, containing the fitted model of a 17 item questionnaire, with 400 observations.
## The model was fitted as linear factor-analytic solution, and it contains 4 doublets
## (items 3-10, 12-13, 14-16 and 11-15) :

SINRELEF.LD(L = FIQ_CONT_L, PSI = FIQ_CONT_VAR, doublet_list = FIQ_CONT_ITEMS_D, cor_doublet =
  FIQ_CONT_VALUES_D, model = "linear", N=400)

## The original subset can be found in:
data(FIQ400)
}
