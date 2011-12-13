\name{rpt.mcmcLMM.adj}
\alias{rpt.mcmcLMM.adj}
\title{
	LMM-based adjusted repeatability estimated using MCMC sampling
}
\description{
	Calculates (adjusted) repeatability from a linear mixed-effects models fitted by MCMC
}
\usage{
	rpt.mcmcLMM.adj(formula, grname, data, CI=0.95, prior=NULL, verbose=FALSE, ...)
}
\arguments{
	\item{formula}{Formula as used e.g. by \link{lmer}. The grouping factor of interest needs to be included as a random effect, e.g. '(1|groups)'. Covariates and additional random effects can be included to estimate adjusted repeatabilities.}
	\item{grname}{A character string or vector of character strings giving the name(s) of the grouping factor(s), for which the repeatability should be estimated. Spelling needs to match the random effect names as given in \code{fromula}.}	
	\item{data}{A dataframe that contains the variables included in the formula argument.}
	\item{CI}{Width of the Bayesian credible interval (defaults to 0.95)}
	\item{prior}{List of prior values passed to the \link{MCMCglmm} function in \pkg{MCMCglmm} (see there for more details). Currently, the default MCMCglmm priors will be used if prior is \code{NULL}.} 
	\item{verbose}{Whether or not \link{MCMCglmm} should print MH diagnostics are printed to screen. Defaults to FALSE.}
	\item{...}{Additonal arguements that are passed on to \link{MCMCglmm} (e.g. length of chain, thinning interval).}
}
\value{
	Returns an object of class rpt that is a a list with the following elements: 
	\item{datatype}{Response distribution (here: "Gaussian").}
	\item{method}{Method used to calculate repeatability (intra-class correlation, ICC) (here: "MCMC").}
	\item{CI}{Width of the Bayesian credibility interval.}
	\item{R}{Point estimate for repeatability (intra-class correlation, ICC), i.e. the mode of the posterior distribution.}
	\item{se}{Standard error (\emph{se}) for repeatability (ICC), i.e. the standard deviation of the posterior distribution. Note that the distribution might not be symmetrical, in which case the se is less informative.}
	\item{CI.R}{Bayesian credibility interval for the repeatability (ICC) based on the posterior distribution of \emph{R}.}
	\item{P}{Significace test, returned as  \code{NA}, since the Bayesian approach conflicts with the null hypothesis testing.}
	\item{R.post}{MCMC samples form the posterior distributions of \emph{R}.} 
}
\details{
	Models are fitted using the \link{MCMCglmm} function in \pkg{MCMCglmm}. 
}
\references{
	Carrasco, J. L. and Jover, L.  (2003). \emph{Estimating the generalized concordance correlation coefficient through variance components}. Biometrics 59: 849-858.
	
	Nakagawa, S. and Schielzeth, H. (2010) \emph{Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists}. Biological Reviews 85: 935-956
}
\author{
	Holger Schielzeth  (holger.schielzeth@ebc.uu.se) & Shinichi Nakagawa (shinichi.nakagawa@otago.ac.nz)
}
\seealso{
	\link{rpt.mcmcLMM}, \link{rpt.remlLMM.adj}, \link{rpt.adj}, \link{print.rpt}
}
\examples{
# repeatability estimation for weight (body mass) - a lower R than the previous one
#data(BodySize)
#(rpt.Weight <- rpt.mcmcLMM.adj(Weight ~ Tarsus + (1|BirdID), "BirdID", data=BodySize))
# A demonstration that the repeatability can be calculated for multiple levels simultanously. In this case the second component is the residuals variance, i.e. the non-repeatable component
#data(BodySize)
#(rpt.Weight <- rpt.mcmcLMM.adj(Weight ~ Tarsus + (1|BirdID), c("BirdID", "units"), data=BodySize))
}
\keyword{models}