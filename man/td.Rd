\name{td}
\alias{td}
\alias{CalcBeta}
\alias{CalcC}
\alias{CalcLogL}
\alias{CalcPowerMatrix}
\alias{CalcQ}
\alias{CalcQ_Lit}
\alias{CalcR}
\alias{CalcRSS}
\alias{CalcS}
\alias{CalcSigma2}
\alias{SubDenton}
\alias{SubRegressionBased}

\title{Temporal Disaggregation of Time Series}
\description{Performs temporal disaggregation and interpolation of low frequency to high frequency time series. \code{td} only works with \code{mts} or \code{ts} objects.}
\usage{
td(formula, conversion = "sum", method = "chow-lin-minrss", 
  to = "quarterly", start = NULL, end = NULL, ...)
}
\arguments{
  \item{formula}{an object of class \code{"\link{formula}"}: a symbolic description of the the temporal disaggregation model. The details of model specification are given under 'Details'.}
  \item{conversion}{type of conversion: \code{"sum"} (default), \code{"average"}, \code{"first"} or \code{"last"}.}
  \item{method}{method of temporal disaggregation: \code{"chow-lin-minrss"} (default), \code{"chow-lin-maxlog"}, \code{"fernandez"}, \code{"litterman"}, \code{"denton-cholette"}, \code{"denton"}, \code{"uniform"} or \code{"ols"}.}
  \item{to}{(high-frequency) destination frequency as a character string (\code{"quarterly"} (default) or \code{"monthly"}) or as a scalar (e.g. \code{2}, \code{4}, \code{12}). Only set if method is equal to  \code{"denton-cholette"}, \code{"denton"} or \code{"uniform"} and no indicators are supplied.}
  \item{start}{(optional) start date. Similar to pre-processing the input series with \code{window}.}
  \item{end}{(optional) end date. Similar to pre-processing the input series with \code{window}.}
  \item{...}{method specific options. See 'Details'.}
}
\details{
\code{td} is used to disaggregate and interpolate a low frequency to a higher frequency time series, while either the sum, the average, the first or the last value of the resulting high-frequency series is consistent with the low frequency series. Disaggregation can be performed with or without the help of one or more indicator series.

The selection of a temporal disaggregation model is similar to the selection of a linear regression model. Thus, \code{td} closely mirrors the working of the \code{\link{lm}} function. The left hand side of the \code{\link{formula}} denotes the low-frequency series, the right hand side the indicators. If no indicator is specified, the right hand side must be set equal to \code{0} (see examples). Unlike \code{lm}, \code{td} handles \code{\link{ts}} and \code{mts} time-series objects, as a typical application involves the use of these objects.

For the generalized least squares (GLS) methods \code{"chow-lin-minrss"}, \code{"chow-lin-maxlog"},  \code{"litterman-minrss"} and \code{"litterman-maxlog"}, an autoregressive parameter \code{rho} is estimated. There are two variants of the \code{chow-lin-minrss} approach in other software packages that lead to different results: Ecotrim by R. Barcellan uses a correlation matrix instead of the variance covariance matrix, the Matlab library by E.M. Quilis multiplies the correlation matrix with 1/(1-rho^2). The two implementations can be replicated by setting \code{vcov} equal to \code{"ecotrim"} (default) or \code{"quilis"}.

The denton methods \code{"denton"} and \code{"denton-cholette"} can be specified with one or without an indicator. The parameter \code{h} can be set equal to \code{0}, \code{1} (default), or \code{2}. Depending on the value, the \code{denton} procedure minimizes the sum of squares of the deviations between the levels (\code{0}), the first differences (\code{1}) or the second differences (\code{2}) of the indicator and the resulting series. Additionally, \code{criterion} can be set equal to \code{"proportional"} (default) or \code{"additive"}, depending on whether the proportional or the absolute deviations should be considered for minimzation. \code{"denton-cholette"} removes the transient movement of the original \code{"denton"} method at the beginning of the resulting series.

\code{"uniform"} is a special case of the \code{"denton"} approach, with \code{h} equals  \code{0} and \code{criterion} equals  \code{"proportional"}. It distributes the residuals uniformly. If no indicator is used, this leads to a step-shaped series.

\code{"ols"} performs an ordinary least squares regression (OLS) and distributes the residuals uniformly. It is useful for comparing the estimators of GLS and OLS regressions.
}

\value{
    \code{td} returns an object of class \code{"td"}.
  
The function \code{predict} computes the interpolated high frequency series. If the high-frequency indicator series are longer than the low-frequency series, the resulting series will be extrapolated. The function \code{coefficients} extracts the coefficients. The function \code{residuals} extracts the low frequency residuals.
  
  An object of class \code{"td"} is a list containing the following components:
  \item{fitted.values}{interpolated (and extrapolated) high frequency series}
  \item{p}{preliminary high frequency series}
  \item{residuals}{low-frequency residuals}
  \item{rho}{autoregressive parameter}
  \item{coefficients}{a named vector of coefficients}
  \item{se}{standard errors of the coefficients}
  \item{s_2}{ML-estimator of the variance of the high-frequency residuals}  
  \item{s_2_gls}{GLS-estimator of the variance of the high-frequency residuals}
  \item{tss}{weighted (low frequency) total sum of squares}
  \item{rss}{weighted (low frequency) residual sum of squares}
  \item{logl}{log-likelihood}
  \item{rank}{number of right hand variables (including intercept)}
  \item{df}{degrees of freedom}
  \item{method}{method of temporal disaggregation}
  \item{call}{function call}
  \item{name}{name of the low frequency variable}
  \item{fr}{the ratio of high- to low-frequency series}
  \item{conversion}{type of conversion}
  \item{actual}{actual values of the low-frequeny series}
  \item{model}{a matrix containing the indicators (and a constant if present)}
  \item{vcov}{type of the variance-covariance approx. in min RSS methods}
  \item{criterion}{minimization criterion in denton methods}
  \item{h}{order of differencing in denton methods}
}

\references{
Chow, G.C. and Lin, A., The Review of Economics and Statistics, Vol. 53, 372--375, 1971

Barcellan, R., ECOTRIM: Software for Temporal Disaggregation, \url{http://circa.europa.eu/Public/irc/dsis/ecotrim/library}

Quilis, E.M. Matlab Temporal Disaggregation Library,
\url{http://www.mathworks.com/matlabcentral/fileexchange/24438-temporal-disaggregation-library}
}


\author{
Christoph Sax, Peter Steiner
}


\seealso{
\code{\link{ta}} for temporal aggregation, the inverse function of \code{td}.

\code{\link[=summary.td]{summary}} is used to obtain and print a summary of the results.
 
\code{\link[=plot.td]{plot}} is used to plot the predicted and actual low frequency series, as well as the residuals.
}
\examples{
data(swisspharma)

# series to be temporally disaggregated
plot(sales.a, main = "Annual Sales, Pharmaceuticals and Chemicals")

# indicator series
plot(exports.q, main = "Quarterly Exports, Pharmaceuticals and Chemicals")

# one indicator, no intercept
mod1 <- td(sales.a ~ 0 + exports.q)
summary(mod1)
plot(mod1)

# disaggregated and true series
ts.plot(cbind(predict(mod1), sales.q), col = c("red", "blue"),
  main = "Quarterly Sales, Pharmaceuticals and Chemicals")
legend("topleft", c("Estimated values", "True values"), 
  lty = 1, col = c("red", "blue"), inset = 0.05)

# temporally aggregated series is equal to the annual value
ta(predict(mod1), conversion = "sum", to = "annual") - sales.a

# several indicators, including intercept
mod2 <- td(sales.a ~ imports.q + exports.q)
summary(mod2)  
plot(predict(mod2), main = "Quarterly Sales, Pharmaceuticals and Chemicals")
  
# no indicator (Denton-Cholette) 
mod3 <- td(sales.a ~ 0, to = "quarterly", method = "denton-cholette", 
  criterion = "proportional", h=2)
plot(predict(mod3), main = "Quarterly Sales, Pharmaceuticals and Chemicals")

# no indicator (uniform)
mod4 <- td(sales.a ~ 0, to = "quarterly", method = "uniform")
plot(predict(mod4), main = "Quarterly Sales, Pharmaceuticals and Chemicals")

# monthly disaggregation
mod.m <- td(sales.a ~ exports.m)
summary(mod.m)
plot(predict(mod.m), main = "Monthly Sales, Pharmaceuticals and Chemicals")

# Quilis vs. Ecotrim
mod2ecotrim <- td(sales.a ~ imports.q + exports.q, vcov = "ecotrim")
mod2quilis <- td(sales.a ~ imports.q + exports.q, vcov = "quilis")
summary(mod2ecotrim)
summary(mod2quilis)
}

\keyword{disaggregation method}
