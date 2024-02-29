\name{SINRELEF.LD-package}
\alias{SINRELEF.LD-package}
\title{
Scores Information, reliability, and relative efficiency in measures that contain locally-dependent items.
}
\description{
\code{SINRELEF.LD} implements an approach aimed at assessing the accuracy and effectiveness of raw scores obtained in scales that contain locally dependent items. The program uses as input the calibration (structural) item estimates obtained from fitting extended unidimensional factor-analytic solutions in which the existing local dependencies are included. Measures of reliability (Omega) and information are proposed at three levels: (a) total score, (b) bivariate-doublet, and (c) item-by-item deletion, and  are compared to those that would be obtained if all the items had been locally independent. All the implemented procedures can be obtained from: (a) linear factor-analytic solutions in which the item scores are treated as approximately continuous, and (b) non-linear solutions in which the item scores are treated as ordered-categorical.
}
\details{
For more information about the methods used, please go to main page \code{\link{SINRELEF.LD}}.
}
\value{
\item{\code{\link{SINRELEF.LD}}}{Assesses the accuracy and effectiveness of raw item and test scores, using the calibration estimates obtained from fitting unidimensional factor-analytic solutions in which item local dependences are included.}
}
\references{

Ferrando, P. J., & Morales-Vives, F. (2023). Is it quality, is it redundancy, or is model inadequacy? Some strategies for judging the appropriateness of high-discrimination items. Anales de Psicologia, 39(3), 517. \doi{10.6018/analesps.535781}

Ferrando, P. J., Navarro-Gonzalez, D., & Lorenzo-Seva, U. (2019). Assessing the Quality and Effectiveness of the Factor Score Estimates in Psychometric Factor-Analytic Applications. Methodology, 15(3), 119-127. \doi{10.1027/1614-2241/a000170}

Raykov, T. (2001). Estimation of congeneric scale reliability using covariance structure analysis with nonlinear constraints. British Journal of Mathematical and Statistical Psychology, 54(2), 315-323. \doi{10.1348/000711001159582}
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
}
\keyword{ package }
