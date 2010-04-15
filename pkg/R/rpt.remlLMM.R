rpt.remlLMM <- function(y, groups, CI=0.95, nboot=1000, npermut=1000) {
	# initial checks
	if(length(y)!= length(groups)) stop("y and group are of unequal length")
	if(class(groups)!="factor") {
		warning("groups will be converted to a factor")
		groups <- factor(groups) 
	}
	if(any(is.na(y))) {
		warning("missing values in y are removed")
		groups <- groups[!is.na(y)]
		groups <- factor(groups)
		y      <- y[!is.na(y)] 
	}
	# preparation
    k <- length(unique(groups))
	N <- length(y)
	# functions: point estimates of R
	R.pe <- function(y, groups) {
		varComps <- nlme::VarCorr(lme(y ~ 1, random = ~ 1|groups ))
		var.a    <- as.numeric(varComps[1,1])
		var.e    <- as.numeric(varComps[2,1])
		R        <- var.a / (var.a + var.e)
		return(R) 
	}
	# point estimation according to model 8 and equation 9
	R <- R.pe(y, groups)
	# confidence interval estimation by parametric bootstrapping
	bootstr <- function(y, groups, k, N, beta0, var.a, var.e) {
		y.boot <- rep(rnorm(1, beta0[1], beta0[2]),N) + rnorm(k, 0, sqrt(var.a))[groups] + rnorm(N, 0, sqrt(var.e))
		R.pe(y.boot, groups) 
	}
	mod      <- lme(y ~ 1, random = ~ 1|groups )
	beta0    <- as.numeric(summary(mod)$tTable[,1:2])
	varComps <- nlme::VarCorr(mod)
	var.a    <- as.numeric(varComps[1,1])
	var.e    <- as.numeric(varComps[2,1])	
	R.boot   <- replicate(nboot, bootstr(y, groups, k, N, beta0, var.a, var.e), simplify=TRUE) 	
	CI.R     <- quantile(R.boot, c((1-CI)/2,1-(1-CI)/2))
	se       <- sd(R.boot)
	# significance test by likelihood-ratio-test
	LR       <- as.numeric(-2*(logLik(lm(y~1))-logLik(mod)))
	P.LRT    <- ifelse(LR<=0, 1, pchisq(LR,1,lower.tail=FALSE)/2)
	# significance test by permutation
	permut <- function(y, groups, N) {
		samp <- sample(1:N, N)
		R.pe(y, groups[samp]) 
	}
	R.permut <- replicate(npermut, permut(y, groups, N), simplify=TRUE)
	P.permut <- sum(R.permut >= R)/npermut
	# return of results
	res  <- list(datatype="Gaussian", method="LMM.REML", CI=CI, R=R, se=se, CI.R=CI.R, P = c(P.LRT=P.LRT, P.permut=P.permut), R.boot=R.boot, R.permut=R.permut )
	class(res) <- "rpt"
	return(res) 
}
