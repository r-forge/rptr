rpt.remlLMM.adj = function(formula, grname, data, CI=0.95, nboot=1000, npermut=1000) {
	mod         <- lmer(formula, data=data)
	# point estimates of R
	R.pe <- function(formula, data, grname) {
		mod.fnc = lmer(formula, data)
		varComps <- lme4::VarCorr(mod.fnc)
		var.a    <- as.numeric(varComps[grname])
		var.p    <- sum(as.numeric(summary(mod.fnc)@REmat[,"Variance"]))
		#var.e    <- as.numeric(attr(varComps, "sc")^2)
		R        <- var.a / var.p
		return(R) 
	}
	R <- R.pe(formula, data, grname)
	names(R) = grname
	# confidence interval estimation by parametric bootstrapping
	bootstr <- function(mod, formula, data, grname) {
		mod.sim <- arm::sim(mod, n.sim=2)   # for some reason it is not possible to set n.sim=1
		y <- mod@X %*% matrix(mod.sim@fixef[1,]) + t(mod@Zt) %*% matrix(unlist(mod.sim@ranef)[seq(1, length(unlist(mod.sim@ranef)), by=2)])
		y <- y + rnorm(length(y), 0, attr(lme4::VarCorr(mod), "sc"))
		data[,names(mod@frame)[1]] = as.vector(y)
		R.pe(formula, data, grname)
	}
	if(any(R==0)) {
		warning("(One of) the point(s) estimate for the repeatability was exactly zero; parametric bootstrapping has been skipped.")
		R.boot = R.pe
		CI.R = rep(NA, length(R.pe))
		names(CI.R) = grname
		se <- rep(NA, length(R.pe))
		names(se) = grname 
	}
	else {
		R.boot   <- replicate(nboot, bootstr(mod, formula, data, grname), simplify=TRUE)
		if(length(grname) == 1) {
			CI.R     <- quantile(R.boot, c((1-CI)/2,1-(1-CI)/2))
			se <- sd(R.boot)
			names(se) = grname }
		else {
			CI.R     <- t(apply(R.boot, 1, function(x) { quantile(x, c((1-CI)/2,1-(1-CI)/2))}))	
			se       <- apply(R.boot, 1, sd)
			rownames(R.boot) = grname
			rownames(CI.R) = grname
			names(se) = grname	}
	}
	# significance test by permutation
	P.permut <- NA
	# significance test by likelihood-ratio-test
	P.LRT    <- NA
	# preparing results
	res  = list(datatype="Gaussian", 
				method="LMM.REML", 
				CI=CI,
				R=R, se=se, CI.R=CI.R, 
				P = c(P.LRT=P.LRT, P.permut=P.permut),
				R.boot=R.boot, R.permut=NA,
				mod = mod	)
	class(res) <- "rpt"
	return(res)
}
