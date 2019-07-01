###############################################################################
# Functions used in the simulation of macroparasite populations. 
# 
###############################################################################

## source file containing functions that are used in parallel implementations 
## and therefore need to be sourced within the parallel ref frame. 
source("ParallelFuncs.R")

## Get parameter value for the name param in the given file. 
## Format of the file:
## paramName (tab) value1 value2 value3...valueN (tab) Comments. 
## returns a vector of strings to be converted if necessary. 
readParam <- function(paramName,fullFilePath)
{
	con <- file(fullFilePath,open="r")
	
	value <- NA
	found = FALSE
	while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0 & !found)
	{
		qq <- strsplit(line,"\t")
		tokens <- qq[[1]]
		if(length(tokens)>0)
		{
			if(tokens[1]==paramName)
			{
				#value <- as.numeric(tokens[2])
				value <- strsplit(tokens[2]," ")
				found=TRUE
			}
		}
	}
	
	close(con)
	if(!found) return(NA) 
	return(value[[1]])
}

## get value of mean worm burden associated with this mean epg and the exponential model. 
expoBurdenFromEPG <- function(testEPG,p)
{
	uniroot(function(m,p1,EPG) epgPerPerson(m,p1)-EPG,interval=c(0,1/(1-p$z)), p1=p,EPG=testEPG)$root
}

############################################################################### 
### functions for power-law dependency (OR OTHER) on worm burden. 
## the relationship between epg/female worm and female worm burden. 
## returns a vector of the log(values) corresponding to the worm burdens found in currentM. 
logPowerLawFecundFunc <- function(currentM, p)
{
	return(log(p$lambdaPL) + log(currentM)*p$exponentPL)	
}

## interpolates an EPG (fertilized and unfertilized) value for testM by interpolation from powerLawEPG list.
## takes the powerLawEPG structure. 
powerLawAllEPGfn <- function(testM,plEPG)
{
	index <- testM%/%plEPG$dM + 1
	offset <- (testM%%plEPG$dM)/plEPG$dM
	return(plEPG$all[index+1]*offset + plEPG$all[index]*(1-offset))
}

## interpolates an EPG (fertilized) value for testM by interpolation from powerLawEPG list.
## takes the powerLawEPG structure. 
powerLawFertileEPGfn <- function(testM,plEPG)
{
	index <- testM%/%plEPG$dM + 1
	offset <- (testM%%plEPG$dM)/plEPG$dM
	return(plEPG$fertile[index+1]*offset + plEPG$fertile[index]*(1-offset))
}

## get value of mean worm burden associated with this mean epg and the power law model.
## the upper limit on worm burden in the function is arbitrary. 
powerLawBurdenFromEPG <- function(testEPG,plEPG)
{
	uniroot(function(m,EPG,plEPG) powerLawAllEPGfn(m,plEPG)-EPG,interval=range(plEPG$interpM),EPG=testEPG,plEPG=plEPG)$root
}


## set up the list structure containing the info for power law dependency of egg output. 
## dM : interval in worm burden between function evaluations. 
## interpMaxM : the max value of worm burden we go up to.
## maxN : max value for the internal integral. 
## Output is the powerLawEPG structure. 
# try: dM=0.1 or 0.05, interpMaxM = 15, maxN=200
powerLawEPGSetUp <- function(dM,interpMaxM,maxN)
{
	## the theory for this is in your notes. 
	interpM <- seq(0,interpMaxM,by=dM)	## Do I actually need this?? 
	interpLength <- length(interpM)
	
	myInts <- 1:maxN
	logInts <- log(myInts)	## to use in the calculation in calculateInterpValues(). 
	logFecDensity <- logPowerLawFecundFunc(1:maxN,params)
	powerLawEPG <- list(interpM=interpM,all=rep(0,interpLength),fertile=rep(0,interpLength),dM=dM) ## interpolation values: default=0. 
	
	powerLawEPG$all[1]=0
	powerLawEPG$fertile[1]=0
	
	## loop through the M values. 
	for(i in 2:interpLength)
	{
		logElements <- logInts + logFecDensity + dnbinom(myInts,mu=powerLawEPG$interpM[i],size=params$k,log=TRUE)
		powerLawEPG$all[i] <- sum(exp(logElements))
		reproFactor <- log( 1 - ((params$k+powerLawEPG$interpM[i])/(params$k+2*powerLawEPG$interpM[i]))^(myInts + params$k) )
		powerLawEPG$fertile[i] <- sum(exp(logElements + reproFactor))
	}
	return(powerLawEPG)
}



## New calculation of R0 with higher resolution too. 
getPsi_bParasites <- function(pCurrent)
{
	## higher resolution. 
	deltaT <- 0.1
	
	muBreaks <- c(0, pCurrent$upperBoundData) 
	
	## inteval-centered ages for the age intervals. 
	modelAges <- seq(min(pCurrent$ageBreaks),max(pCurrent$ageBreaks)-deltaT,by=deltaT) + 0.5*deltaT ## each age group in the actual model is annual characterised by mid-value (0.5 is [0,1), etc.)
	
	## hostMu for the new age intervals. 
	hostMuGroupIndex <- cut(modelAges,breaks=muBreaks,labels=1:length(pCurrent$hostMuData)) 
	hostMu <- pCurrent$hostMuData[hostMuGroupIndex] 
	
	meanDeaths <- hostMu*deltaT
	hostSurvivalCurve <- exp(-cumsum(meanDeaths))
	
	hostSurvivalTotal <- sum(hostSurvivalCurve[1:length(modelAges)])*deltaT  ## This should be the integral under the survival curve UP TO THE TOP AGE LIMIT. 
	
	## calculate the cumulative sum of host and worm death rates from which to calculate worm survival. 
	intMeanWormDeathEvents <- cumsum(hostMu+pCurrent$sigma)*deltaT
	
	## need rho and beta at this age resolution as well. 
	ageIndices <- 1:(length(pCurrent$ageBreaks)-1) 
	modelAgeGroupCatIndex <- cut(modelAges,breaks=pCurrent$ageBreaks,labels=ageIndices) 
	betaAge <- pCurrent$betaGroup[modelAgeGroupCatIndex]
	rhoAge <- pCurrent$rhoGroup[modelAgeGroupCatIndex] 
	
	
	## calculate the infectiousness bit first. 
	nMax <- length(hostMu)
	K <- rep(0,nMax)

	for(i in 1:nMax)
	{
		a <- i:nMax
		currentMeanDeath <- intMeanWormDeathEvents[a] - intMeanWormDeathEvents[i]
		K[i] <- deltaT*sum(rhoAge[a]*exp(-currentMeanDeath))
	}
	
	summation <- sum(betaAge*hostSurvivalCurve*K)*deltaT
	
	psi <- pCurrent$R0*hostSurvivalTotal*pCurrent$lDecayRate/(pCurrent$lambda*pCurrent$z*summation)
	
	return(psi)
}


## mean, weighted by population demography.
# x - the quantity to be averaged by age category. 
# survival - the survival curve. 
# ageCats - a vector associating the ith annual age category with the label of each age category (e.g. 1 - infants, 2 - preSAC, etc)
getMean <- function(x,survival,ageCats)
{
	## calculate means in age categories...
	meanSurvival <- x*survival
	meanSurvivalByCat <- tapply(meanSurvival,ageCats,sum)
	normByCat <- tapply(survival,ageCats,sum)
	return(meanSurvivalByCat/normByCat)
}


## For the age structured model used in AgeStructSeparatrix.R and AgeStructureModel3.R
## Calculate and return the correct R0 for the given parameter values.  
getR0 <- function(pCurrent)
{
	return(NA) ## this needs to rewritten to reflect changes to the calculation of R0. 
}

## For the age structured model used in AgeStructSeparatrix.R and AgeStructureModel3.R
## a derivative function to return the value of the derivatives at the given time and state. 
ageStructDerivs <- function(t,state,pCurrent)  #  used to have this: ,derivativeFuncName
{
	M <- state[1:pCurrent$nAgeGroups] 
	L <- state[pCurrent$nAgeGroups+1] 
	
	## change in M... 
	gain <- pCurrent$betaAge*L
	loss <- -(pCurrent$sigma + 1/pCurrent$wAgeGroup)*M 		## Removed the term pCurrent$hostMu from the loss rates, since it's already included in the survival function.
	aging <- c(0,M[1:(pCurrent$nAgeGroups-1)])/pCurrent$wAgeGroup
	dM <- loss + aging + gain
	
	# Change in L... 
	## apply the epg function to all elements of M. 
	#totalInput <- pCurrent$psi*sapply(X=M,FUN=pCurrent$reproFunc,p=pCurrent)*pCurrent$hostSurvivalCurve*pCurrent$rhoAge*pCurrent$wAgeGroup/pCurrent$hostSurvivalTotal  ## output of epg * fraction of population * strength of contribution. 
	#totalInput <- pCurrent$psi*pCurrent$reproFunc(M,pCurrent)*pCurrent$hostSurvivalCurve*pCurrent$rhoAge*pCurrent$wAgeGroup/pCurrent$hostSurvivalTotal
	totalInput <- pCurrent$reproFunc(M,pCurrent)*pCurrent$ageStructPreMult2
	dL <- pCurrent$psi*sum(totalInput) - pCurrent$lDecayRate*L
	return(list(c(dM,dL)))
}


## For the age structured model used in AgeStructSeparatrix.R and AgeStructureModel3.R
## a derivative function to return the value of the derivatives at the given time and state, assuming L is at equilibrium. 
ageStructDerivsNoL <- function(t,state,pCurrent)
{
	M <- state[1:pCurrent$nAgeGroups] 
	#totalInput <- pCurrent$psi*sapply(X=M,FUN=pCurrent$reproFunc,p=pCurrent)*pCurrent$hostSurvivalCurve*pCurrent$rhoAge*pCurrent$wAgeGroup/pCurrent$hostSurvivalTotal  ## output of epg * fraction of population * strength of contribution. 
	#totalInput <- pCurrent$psi*pCurrent$reproFunc(M,pCurrent)*pCurrent$hostSurvivalCurve*pCurrent$rhoAge*pCurrent$wAgeGroup/pCurrent$hostSurvivalTotal  ## Didn't need the sapply because repro functions are vector!!
	totalInput <- pCurrent$reproFunc(M,pCurrent)*pCurrent$ageStructPreMult2
	L <- pCurrent$psi*sum(totalInput)/pCurrent$lDecayRate
	
	## change in M... 
	gain <- pCurrent$betaAge*L
	loss <- -(pCurrent$sigma + 1/pCurrent$wAgeGroup)*M			## Removed the term pCurrent$hostMu from the loss rates, since it's already included in the survival function. 
	aging <- c(0,M[1:(pCurrent$nAgeGroups-1)])/pCurrent$wAgeGroup
	dM <- loss + aging + gain
	
	# Change in L... 
	## apply the epg function to all elements of M. 
	#totalInput <- pCurrent$psi*sapply(X=M,FUN=derivativeFuncName,p=pCurrent)*pCurrent$hostSurvivalCurve*pCurrent$rhoAge*pCurrent$wAgeGroup/pCurrent$hostSurvivalTotal  ## output of epg * fraction of population * strength of contribution. 
	#dL <- sum(totalInput) - pCurrent$lDecayRate*L
	return(list(c(dM)))
}


## Takes a matrix of mean worm burdens with time rows and age columns and returns a matrix with fraction in high worm burden category. 
getFracHighBurdenByAge <- function(MValues)
{
	getHighFrac <- function(x) 1 - pnbinom(params$highBurdenAge,size=params$k,mu=x)
	heavyBurdenFraction <- t(apply(MValues,MARGIN=1,getHighFrac))
	return(heavyBurdenFraction)
}

########################################################################################################################
### High EPG burden fraction...

## calculate the the probability of mean worm burden generating an EPG greater than the threshold. 
epgThreshold3 <- function(currentMu)
{
	## where's the top end for this mu and k? 
	top <- qnbinom(0.995,mu=currentMu,size=params$k)
	nf <- 1:top
	p_nf <- dnbinom(nf,mu=currentMu,size=params$k)
	meanEPG <- params$lambda*nf*params$z^nf
	pGreaterEBvsNf <- 1 - pnbinom(params$highEPGThreshold,mu=meanEPG,size=params$k_epg)
	pGreaterEB <- sum(p_nf*pGreaterEBvsNf)
	return(pGreaterEB)
}

## Takes a vector of M values and returns an array of P(E>Eb)
epgThreshold2 <- function(x)
{
	probsGraterEb <- sapply(x,FUN=epgThreshold3)
	return(probsGraterEb)
}

## Takes a matrix of mean worm burdens with time rows and age columns and returns a matrix with fraction in high EPG category. 
epgThreshold <- function(MValues)
{
	heavyEPGFraction <- t(apply(MValues,MARGIN=1,epgThreshold2))
	return(heavyEPGFraction)
}

inputParameters <- function(paramFilePath,demogName="Default")
{
	params <- list()
	## Epidemiological parameters.. 
	params$lambda <- as.numeric(readParam("lambda",paramFilePath)) # eggs / gram [See below]
	params$z <- exp(-as.numeric(readParam("gamma",paramFilePath))) 		## JET_Nov13
	params$R0 <- as.numeric(readParam("R0",paramFilePath))
	params$lDecayRate <- as.numeric(readParam("ReservoirDecayRate",paramFilePath))
	params$sigma <- as.numeric(readParam("sigma",paramFilePath))
	params$k <- as.numeric(readParam("k",paramFilePath)) ## population k.  
	params$psi <- 1 ## dummy value prior to R0 calculation. 
	
	
	## 'external' age structure.  
	#params$ageBreaks <- c(0,10,17,31,maxIndAge)  ## integer values for the age categories of interest. 
	params$ageBreaks <- as.numeric(readParam("contactAgeBreaks",paramFilePath))  ## integer values for the age categories of interest. 
	params$betaGroup <- as.numeric(readParam("betaValues",paramFilePath))
	params$rhoGroup <- as.numeric(readParam("rhoValues",paramFilePath))
	
	## treatment breaks: these aren't the same as contact structure defined above. 
	params$treatmentBreaks <- as.numeric(readParam("treatmentBreaks",paramFilePath))
	params$coverage <- as.numeric(readParam("coverage",paramFilePath))
	params$drugEff <- as.numeric(readParam("drugEff",paramFilePath))
	params$treatStart <- as.numeric(readParam("treatStart",paramFilePath))
	params$treatmentRounds <- as.numeric(readParam("nRounds",paramFilePath))
	params$interval <- as.numeric(readParam("treatInterval",paramFilePath))  ## interval between treatments. 
	
	## assign reproduction function objected to list of params. 
	params$reproFuncName <- readParam("reproFuncName",paramFilePath) ## for the record...
	params$reproFunc <- match.fun(params$reproFuncName) ## carry the reproduction function with the parameters. 
	
	## high burden levels with age. 
	params$highBurdenBreaks <- as.numeric(readParam("highBurdenBreaks",paramFilePath))
	params$highBurdenValues <- as.numeric(readParam("highBurdenValues",paramFilePath))
	
	params$highEPGThreshold <- as.numeric(readParam("highEPGThreshold",paramFilePath))
	params$k_epg <- as.numeric(readParam("k_epg",paramFilePath))
	
	##################################################################################################
	## construct the path for the demography file from the param file path. 
	g <- gregexpr("[/\\]",paramFilePath)  ## Using regular expressions to capture the position of slashes in both directions. 
	stub <- substr(paramFilePath,1,max(g[[1]]))
	demogPath <- paste0(stub,"Demographies.txt")
	
	## record the demography used. 
	params$demogType <- demogName 
	
	## construct the parameter names...
	muName <- paste0(demogName,"_hostMuData")
	boundsName <- paste0(demogName,"_upperBoundData")
	
	
	## new code based on host death rate data from the Global Health Observatory Data Repository. 
	params$hostMuData <- as.numeric(readParam(muName,demogPath)) # c(0.050, 0.006, 0.003, 0.003, 0.003, 0.004, 0.005, 0.007, 0.008, 0.008, 0.009, 0.010, 0.013, 0.019, 0.029, 0.048, 0.077, 0.137, 0.214, 0.316, 0.445, 0.592)
	params$upperBoundData <- as.numeric(readParam(boundsName,demogPath)) # c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110)

	return(params)
}

## calculate the things we need from the raw parameter values. 
configureHostsParasites <- function(params)
{
	## if it's monogamous, run the monogamous reproduction configuration function.
	if(params$reproFuncName=="epgMonog")
	{
		params$monogParams <- monogFertilityConfig(params,N=30)	
	}
	
	## 'internal' age structure. 
	ageClassWidth <- 1
	params$modelAges <- seq(0,max(params$ageBreaks)-ageClassWidth,by=ageClassWidth) + 0.5*ageClassWidth ## each age group in the actual model is annual characterised by mid-value (0.5 is [0,1), etc.)
	params$nAgeGroups <- length(params$modelAges)
	params$wAgeGroup <- ageClassWidth ##params$modelAges[2] - params$modelAges[1]  
	
	## reconcile internal and external. 
	ageIndices <- 1:(length(params$ageBreaks)-1) 
	params$modelAgeGroupCatIndex <- cut(params$modelAges,breaks=params$ageBreaks,labels=ageIndices) 
	params$betaAge <- params$betaGroup[params$modelAgeGroupCatIndex]
	params$rhoAge <- params$rhoGroup[params$modelAgeGroupCatIndex] 
	
	#treatmentLevels <- params$coverage*params$drugEff
	#treatmentIndices <- 1:(length(params$treatmentBreaks)-1)
	#params$modelAgeTreatmentGroupIndex <- cut(params$modelAges, breaks=params$treatmentBreaks, labels=treatmentIndices)
	#params$treatmentAge <- 1 - treatmentLevels[params$modelAgeTreatmentGroupIndex]   ### [Doesn't take into account efficacy] The fraction remaining after treatment. Hence the fractional reduction. 
	
	## high burden levels with age. 
	burdenIndices <- 1:(length(params$highBurdenBreaks)-1)
	params$modelAgeBurdenGroupIndex <- cut(params$modelAges,breaks=params$highBurdenBreaks,labels=burdenIndices)
	params$highBurdenAge <- params$highBurdenValues[params$modelAgeBurdenGroupIndex]
	
	muBreaks <- c(0, params$upperBoundData) 
	
	hostMuGroupIndex <- cut(params$modelAges,breaks=muBreaks,labels=1:length(params$hostMuData)) 
	params$hostMu <- params$hostMuData[hostMuGroupIndex] 
	
	meanDeaths <- params$hostMu*ageClassWidth
	params$hostSurvivalCurve <- exp(-cumsum(meanDeaths))
	
	params$hostSurvivalTotal <- sum(params$hostSurvivalCurve)*params$wAgeGroup
	
	## make a more trapezium rule-type integral.  
	tempLength <- length(params$hostSurvivalCurve)
	tempHostSCurve <- (c(1,params$hostSurvivalCurve[1:(tempLength-1)]) + params$hostSurvivalCurve)/2
	tempHostSTotal <- sum(tempHostSCurve)*params$wAgeGroup
	params$ageStructPreMult2 <- tempHostSCurve*params$rhoAge*params$wAgeGroup/tempHostSTotal
	
	## adjust psi to R0. 
	params$psi <- getPsiParasites(params)
	
	return(params)
}