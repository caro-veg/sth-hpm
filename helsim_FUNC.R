# Function file for helsim_RUN.R

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

readParams <- function (fileName,demogName="Default") 
{
	numReps <- as.integer(readParam("repNum",fileName)) # Number of repetitions
	maxTime <- as.numeric(readParam("nYears",fileName)) # Maximum number of years to run
	nYearsPostTreat <- as.numeric(readParam("nYearsPostTreat",fileName)) # The number of years to continue running after treatment (will override nYears if specified)
	N <- as.integer(readParam("nHosts",fileName)) # Host population size
	R0 <- as.numeric(readParam("R0",fileName)) # Basic reproductive number
	lambda <- as.numeric(readParam("lambda",fileName)) # Eggs per gram
	gamma <- as.numeric(readParam("gamma",fileName)) # Exponential density dependence of parasite adult stage
	k <-as.numeric(readParam("k",fileName)) # Shape parameter of assumed negative binomial distribution of worms amongst host
	sigma <- as.numeric(readParam("sigma",fileName)) # Worm death rate
	LDecayRate <- as.numeric(readParam("ReservoirDecayRate",fileName)) # Decay rate of the infectious material in the environment (reservoir decay rate)
	ContactAGB <- as.numeric(readParam("contactAgeBreaks",fileName)) # Contact age group breaks
	beta <- as.numeric(readParam("betaValues",fileName)) # Contact rates
	rho <- as.numeric(readParam("rhoValues",fileName))	# rho, contribution values. 
	CAG <- as.numeric(readParam("coverage",fileName)) # Coverage of incoming migrants
	DrugEfficacy <- as.numeric(readParam("drugEff",fileName)) # Efficacy of drugs 
	outputFrequency <- as.numeric(readParam("outputFrequency",fileName))
	outputOffset <- as.numeric(readParam("outputOffset",fileName))
	highBurdenBreaks <- as.numeric(readParam("highBurdenBreaks",fileName))
	highBurdenValues <- as.numeric(readParam("highBurdenValues",fileName))
	emigrationTimes <- as.numeric(readParam("emigrationTimes", fileName))
	immigrationTimes <- as.numeric(readParam("immigrationTimes", fileName))
	numberPeopleToMigrate <- as.numeric(readParam("numberPeopleToMigrate", fileName))
	migrantAgeGroup <- as.numeric(readParam("migrantAgeGroup", fileName))
	
	
	## sort out times for treatments. 
	treatmentStart <- as.numeric(readParam("treatStart",fileName))
	nRounds <- as.numeric(readParam("nRounds",fileName)) 
	treatInterval <- as.numeric(readParam("treatInterval",fileName)) 
	chemoTimings <- seq(from=treatmentStart,by=treatInterval,length.out=nRounds)

	## simulation-specific parameters
	nNodes <- as.numeric(readParam("nNodes", fileName))
	maxStep <- as.numeric(readParam("maxStep", fileName))
	seed <- as.numeric(readParam("seed", fileName))

	
	pars <- list (	numReps=numReps, 
					maxTime=maxTime, 
					nYearsPostTreat = nYearsPostTreat,
					N=N,
					R0=R0,
					lambda=lambda,
					gamma=gamma,
					k=k,
					sigma=sigma,
					LDecayRate=LDecayRate,
					DrugEfficacy=DrugEfficacy,
					chemoTimings=chemoTimings,
					contactAgeBreaks=ContactAGB,
					nRounds=nRounds,
					contactRates=beta,
					rho=rho,
					coverage=CAG,
					treatInterval=treatInterval,
					treatmentStart=treatmentStart,
					outputFrequency=outputFrequency,
					outputOffset=outputOffset,
					highBurdenBreaks=highBurdenBreaks,
					highBurdenValues=highBurdenValues,
					emigrationTimes=emigrationTimes,
					immigrationTimes=immigrationTimes,
					numberPeopleToMigrate=numberPeopleToMigrate,
					migrantAgeGroup=migrantAgeGroup,
					nNodes=nNodes,
					maxStep=maxStep,
					seed=seed)

	##################################################################################################
	## construct the path for the demogrphy file from the param file path. 
	g <- gregexpr("[/\\]",fileName)  ## Using regular expressions to capture the position of slashes in both directions. 
	stub <- substr(fileName,1,max(g[[1]]))
	demogPath <- paste0(stub,"Demographies.txt")
	
	## record the demography used. 
	pars$demogType <- demogName 
	
	## construct the parameter names...
	muName <- paste0(demogName,"_hostMuData")
	boundsName <- paste0(demogName,"_upperBoundData")
			
	#### host survival curve.
	## Read in the data. 
	pars$hostMuData <- as.numeric(readParam(muName,demogPath)) 
	pars$muBreaks <- c(0,as.numeric(readParam(boundsName,demogPath)))
	
	## turn SR on or off. 
	pars$SR <- TRUE
	if(readParam("StochSR",fileName)=="FALSE")
	{
		pars$SR <- FALSE
	}

	## read in function name for deterministic repro function. Needed for calculation of equilibrium. 
	pars$reproFuncName <- readParam("reproFuncName", fileName)
	
	# Parasite survival curve stuff
	pars$z <- exp(-pars$gamma) # Fecundity parameter (z)
	pars$psi <- 1.0 # Dummy value prior to R0 calculation
		
	return(pars)
}



configure <- function(pars)
{
  	## this was in the read in function...
  	pars$z <- exp(-pars$gamma) # Fecundity parameter (z)
  
  	## look up the name of the reproduction function. 
  	pars$reproFunc <- match.fun(params$reproFuncName) ## carry the reproduction function with the parameters. 
  

	dT <- 0.1 ## level of descretization for the drawing of lifespans. 
	pars$maxHostAge <- min(max(pars$muBreaks),max(pars$contactAgeBreaks))
	
	#pars$muAges <- seq(0,pars$maxHostAge-dT, by=dT)+0.5*dT
	pars$muAges <- seq(0,max(pars$muBreaks)-dT, by=dT) + 0.5*dT   ### the whole range of ages. Concatenate later. 
	muIndices <- cut(pars$muAges,breaks=pars$muBreaks,labels=1:length(pars$hostMuData)) 
	pars$hostMu <- pars$hostMuData[muIndices] 
	
	## first element gives probability of surviving to end of first year. 
	pars$hostSurvivalCurve <- exp(-cumsum(pars$hostMu)*dT)
	
	## the index for the last age group before the cutoff in this descretization. 
	maxAgeIndex <- which(pars$muAges > pars$maxHostAge)[1] - 1
	
	## the cumulative probability of dying in the ith year.  
	fullHostAgeCumulDistr <- cumsum(dT*pars$hostMu*c(1,pars$hostSurvivalCurve[1:(length(pars$hostSurvivalCurve)-1)]))  ## prob of dying in each age group. 
	#pars$hostAgeCumulDistr <- cumsum(dT*pars$hostMu*c(1,hostSurvivalCurve[1:(length(hostSurvivalCurve)-1)]))
	pars$hostAgeCumulDistr <- c(fullHostAgeCumulDistr[1:(maxAgeIndex-1)],1) ## cumulative probability with cutoff at max age. 	

	# Contact stuff
	pars$contactAgeGroupBreaks <- c(pars$contactAgeBreaks[1:(length(pars$contactAgeBreaks)-1)], pars$maxHostAge)  ##  + dT used to have this to avoid going outside bounds of breaks. Not necessary, I think.    
	
	# Timings stuff
	if (!is.na(pars$nYearsPostTreat)){
	  pars$maxTime <- pars$treatmentStart + (pars$nRounds * pars$treatInterval) + pars$nYearsPostTreat # if nYearsPostTreat is specified, update maxTime to the year required
	}
	if(!is.na(pars$outputFrequency)){
		timings = seq(pars$outputOffset, pars$maxTime, 1/pars$outputFrequency) # if outputFrequency is specified, update the output times accordingly
		timings[timings < 0] = 0 # if there are any times smaller than zero (eg from a negative offset) then reset those to zero
		pars$outTimings = timings
	}  
	if(pars$outTimings[length(pars$outTimings)] != pars$maxTime) pars$outTimings = c(pars$outTimings, pars$maxTime) # always output at the end of simulation
	
  return (pars)
}



# Sets up the simulation to initial conditions based on analytical equilibria
# simplified version for migration model:
# takes parameters and L_equi (equilibrium environmental reservoir state)
# only data stored: 
# list 1: total worm counts, female worm counts, flat to indicate if people are in home village or away
# list 2: birth dates, death dates
# si
# freeLiving: equilibrium environmental reservoir state
# contactAgeGroupIndices

setupSD <-  function(pars) 
{
    	si <- rgamma(pars$N,scale=1/pars$k,shape=pars$k)
  
	# Distribute birth dates such that at time zero have the right distribution of ages
	# Give a sample of ages at the start of the simulation
	# For an exponential distribution, the survival function S(t)=exp(-t/mu)
	# So to sample from this distribution, generate a random number on 0 1 and then invert this function
	
	lifeSpans <- getLifeSpans(pars$N, pars)
	trialBirthDates <- -lifeSpans*runif(pars$N)
	trialDeathDates <- trialBirthDates + lifeSpans 
	
	## Equilibrate the population over 1000 years, in the absence of understanding how to generate it in the first place. 
	communityBurnIn <- 1000
	while(min(trialDeathDates)<communityBurnIn)
	{
		earlyDeath <- which(trialDeathDates < communityBurnIn)
		trialBirthDates[earlyDeath] <- trialDeathDates[earlyDeath]
		trialDeathDates[earlyDeath] <- trialDeathDates[earlyDeath] + getLifeSpans(length(earlyDeath),pars) 
	}
	
	demography <- data.frame(birthDate=trialBirthDates-communityBurnIn, deathDate=trialDeathDates-communityBurnIn)
  
 	# Contactparam index for each host
	contactAgeGroupIndices <- cut(-demography$birthDate, pars$contactAgeGroupBreaks, 1:(length(pars$contactAgeGroupBreaks)-1))
	
	
	### calculate the IC worm burdens here... 
	meanBurdenIndex <- cut(-demography$birthDate,breaks=c(0,pars$equiData$ageValues),labels=1:length(pars$equiData$ageValues))
	means <- si*pars$equiData$hatProfile[meanBurdenIndex]*2  ## factor of 2 because the profile is for female worms in the det system. 
	wTotal <- rpois(pars$N,lambda=means) 
	worms <- data.frame(total=wTotal,female=rbinom(pars$N,size=wTotal,prob=0.5))
	worms$away <- rep(FALSE, pars$N)	# flag indicating if individual i is in home village or away
	stableFreeLiving <- pars$equiData$L_stable*2  ## 2 because all eggs not just female ones
	
	
 	SD <- list (ID=1:pars$N, 
			si=si, 
			worms=worms, 
			freeLiving=stableFreeLiving,
               	demography=demography,
			contactAgeGroupIndices=contactAgeGroupIndices
			)
	
	return(SD)
}




# Output the simulation
outRes <- function (f, r, t, SD) 
{
	cat (file=f, r, t, mean(SD$worms$female),SD$worms$female,SD$freeLiving, "\n", sep="\t") # Adults only
}

# Output the results from multiple runs
outResMulti <- function(f, SD) 
{
	cat (file=f, mean(SD$worms$female),"\t") # Adults only
}

# Calculate the event rates


# Calculate the event rates
calcRatesMig <- function(pars, SD, sourceReservoir, t=0) 
{
	# Only two events: worm death and new worms
	# Worm total death rate
	deathRate <- pars$sigma*sum(SD$worms$total)

	# freeLiving - vector
	currentFL <- c(SD$freeLiving, sourceReservoir)
	hostInfRates <- currentFL[1 + as.numeric(SD$worms$away)]*SD$si*pars$contactRates[SD$contactAgeGroupIndices]
		
	rates <- c (hostInfRates,deathRate)
	return(rates)
}



# Enact an event
doEvent <- function(rates, pars, SD, t=0) 
{
	# Determine which event
	# If it's 1 to N, it's a new worm otherwise, it's a worm death
	event <- which((runif(1)*sum(rates))<cumsum(rates))[1]
  
	if(event==length(rates))
	{
		# It's a worm death... 
		deathIndex <- which((runif(1)*sum(SD$worms$total))<cumsum(SD$worms$total))[1]
		# Is it female? 
		if(runif(1)<SD$worms$female[deathIndex]/SD$worms$total[deathIndex]) 
		    SD$worms$female[deathIndex] <- SD$worms$female[deathIndex] - 1
		SD$worms$total[deathIndex] <- SD$worms$total[deathIndex] - 1
	} 
	else
	{
		# It's a new worm...
		SD$worms$total[event] <- SD$worms$total[event] + 1
		if(runif(1)<0.5) 
			  SD$worms$female[event] <- SD$worms$female[event] + 1
	}
	return(SD)  
} 


# Run processes that need to occur regularly, i.e reincarnating whichever hosts have died recently and
# updating the free living worm population
doRegular <- function(pars, SD, ts)
{
  SD <- doDeath(pars, SD, t)
  SD <-  doFreeLive(pars, SD, ts)
  return (SD)
}


# Update the freeliving populations deterministically, use this function in migration simulations
# when individuals might be "away" and not contribute to the environmental reservoir
doFreeLiveMigration <- function(pars, SD, ts) 
{
  
  # Polygamous reproduction - female worms produce fertilised eggs only if there's at least one male worm around
  if(pars$reproFuncName == "epgFertility")
  {
    	noMales <- SD$worms$total==SD$worms$female
    	productivefemaleworms <- SD$worms$female
    	if(pars$SR)
    	{
      	productivefemaleworms[noMales] <- 0		## Commenting out this line would remove sexual reproduction, I think. 
    	}    
  }


  
  # Monogamous reproduction - only pairs of worms produce eggs
  else if(pars$reproFuncName == "epgMonog")
  {
    	male <- SD$worms$total - SD$worms$female
    	productivefemaleworms <- pmin(male, SD$worms$female)    
  }

  

  eggOutputPerHost <- pars$lambda*productivefemaleworms*exp(-productivefemaleworms*pars$gamma)
  eggOutputPerHost <- eggOutputPerHost*as.numeric(1 - SD$worms$away)
  ## Factor of 2 because psi is rate of prod of female eggs. We want total fertilized eggs for total worms.		
  eggsProdRate <- 2*pars$psi*sum(eggOutputPerHost*pars$rho[SD$contactAgeGroupIndices])/pars$N		
  

  ## dL/dt = K - mu*L has soln: L(0)exp(-mu t) + K*(1-exp(-mu t))/mu.  This is exact if rate of egg production is constant in the timestep.
  expo <- exp(-pars$LDecayRate*ts)
  SD$freeLiving <- SD$freeLiving*expo + eggsProdRate*(1-expo)/pars$LDecayRate
  return(SD)
}





# Note that despite the name, this is run periodically to check whether there have been any deaths but also to update 
# the age group the hosts fall into. 

# version without treatmentAgeGroupIndices for BLANK simulation
# version with treated and untreated village would require treatmentAgeGroupIndices
doDeath <- function(pars, SD, t)  ## , logFile
{
	# Identify the indices of the dead 
	theDead <- which(SD$demography$deathDate < t) 
	if(length(theDead!=0))
	{ # note that prior to fix, if there were no dead, this would return to the caller and not update age categories
		
		# Set the birth date to now
		SD$demography$birthDate[theDead] <- t - 0.001   ## put the birth slightly in the past to ensure age is just positive for categorization. 
		
		# Calculate death dates
		newDeathDays <- t + getLifeSpans(length(theDead),pars)
		SD$demography$deathDate[theDead] <- newDeathDays
		
		# They also need new force of infections (FOIs)
		SD$si[theDead] <- rgamma(length(theDead),scale=1/pars$k,shape=pars$k)
		
		# Kill all their worms 
		SD$worms$total[theDead] <- 0
		SD$worms$female[theDead] <- 0
				
	}
	
	
	# Redo the contact age categories
	SD$contactAgeGroupIndices <- cut(t-SD$demography$birthDate,pars$contactAgeGroupBreaks,1:(length(pars$contactAgeGroupBreaks)-1))

	return(SD)
}



# version with treatmentAgeGroupIndices for simulation with treated and untreated village
doDeathWithTreatAgeGroups <- function(pars, SD, t)  
{
	# Identify the indices of the dead 
	theDead <- which(SD$demography$deathDate < t) 
	if(length(theDead!=0))
	{ # note that prior to fix, if there were no dead, this would return to the caller and not update age categories
		
		# Set the birth date to now
		SD$demography$birthDate[theDead] <- t - 0.001   ## put the birth slightly in the past to ensure age is just positive for categorization. 
		
		# Calculate death dates
		newDeathDays <- t + getLifeSpans(length(theDead),pars)
		SD$demography$deathDate[theDead] <- newDeathDays
		
		# They also need new force of infections (FOIs)
		SD$si[theDead] <- rgamma(length(theDead),scale=1/pars$k,shape=pars$k)
		
		# Kill all their worms 
		SD$worms$total[theDead] <- 0
		SD$worms$female[theDead] <- 0
		
	}
	
	
	# Redo the contact age categories
	SD$contactAgeGroupIndices <- cut(t-SD$demography$birthDate,pars$contactAgeGroupBreaks,1:(length(pars$contactAgeGroupBreaks)-1))
	
	# Work out which treatment age category each host falls into  
	treatmentAgeGroupIndices_new <- cut(t-SD$demography$birthDate,pars$treatmentAgeGroupBreaks,1:(length(pars$treatmentAgeGroupBreaks)-1))	
	# Update the treatment age categories
	SD$treatmentAgeGroupIndices <- treatmentAgeGroupIndices_new
	
	return(SD)
}





# Chemotherapy
doChemo <- function (pars, SD, t) 
{  
  # Decide which individuals are treated - treatment is random
  adherenceFactors <- runif(pars$N)
  attendance <- getAttendance(adherenceFactors, pars$coverage, SD$treatmentAgeGroupIndices) 
  
  # How many worms die? 
  maleWorms <- SD$worms$total[attendance] - SD$worms$female[attendance]
  maleToDie <- rbinom(sum(attendance), maleWorms, pars$DrugEfficacy) 
  femaleToDie <- rbinom(sum(attendance), SD$worms$female[attendance], pars$DrugEfficacy) 
  
  SD$worms$total[attendance] <- SD$worms$total[attendance] - maleToDie - femaleToDie
  SD$worms$female[attendance] <- SD$worms$female[attendance] - femaleToDie
  
  return(SD)  
} 



# Return the value of psi, given the parameters and R0 
getPsiOLD <- function(pars)
{
	qVector <- pars$hostMu + pars$sigma + 1  # The 1 is for movement from one class to the next
	preSurvivalLog <- c(0,-log(1+pars$hostMu+pars$sigma)) 
	wormSurvivalLog <- cumsum(preSurvivalLog) 

	# Construct the matrix and sum it
	summationMatrix <- matrix(0,pars$NAG,pars$NAG)
	ageIndices <- cut(1:pars$NAG-0.1,pars$contactAgeGroupBreaks,1:(length(pars$contactAgeGroupBreaks)-1))
	betaAge <- pars$contactRates[ageIndices]
	for(i in 1:pars$NAG)
	{
		for(j in i:pars$NAG)
		{
			summationMatrix[i,j] <- betaAge[i]*pars$hostSurvivalCurve[j]*exp(wormSurvivalLog[j] - wormSurvivalLog[i])/qVector[j]
		}
	}
	
	summation <- pars$lambda*pars$z*sum(summationMatrix)/(pars$hostSurvivalTotal*pars$LDecayRate) 
	
	summation <- summation * 0.5   # Only half the worms are female and the model is in terms of female worms
	
	return(pars$R0/summation)
}

## Calculate and return the correct psi to give the required R0. 
## New calculation of R0 with higher resolution too. 
getPsi_b <- function(pCurrent)
{
	## higher resolution. 
	deltaT <- 0.1
	
	## inteval-centered ages for the age intervals. 
	## 0 to maxHostAge, mid points. 
	modelAges <- seq(0,pCurrent$maxHostAge-deltaT,by=deltaT) + 0.5*deltaT ## each age group in the actual model is annual characterised by mid-value (0.5 is [0,1), etc.)
	
	## hostMu for the new age intervals. 
	hostMuGroupIndex <- cut(modelAges,breaks=pCurrent$muBreaks,labels=1:length(pCurrent$hostMuData)) 
	hostMu <- pCurrent$hostMuData[hostMuGroupIndex] 
	
	meanDeaths <- hostMu*deltaT
	hostSurvivalCurve <- exp(-cumsum(meanDeaths))
	
	hostSurvivalTotal <- sum(hostSurvivalCurve[1:length(modelAges)])*deltaT  ## This should be the integral under the survival curve UP TO THE TOP AGE LIMIT. 
	
	## calculate the cumulative sum of host and worm death rates from which to calculate worm survival. 
	intMeanWormDeathEvents <- cumsum(hostMu+pCurrent$sigma)*deltaT
	
	## need rho and beta at this age resolution as well. 
	ageIndices <- 1:(length(pCurrent$contactAgeGroupBreaks)-1) 
	modelAgeGroupCatIndex <- cut(modelAges,breaks=pCurrent$contactAgeGroupBreaks,labels=ageIndices) 
	betaAge <- pCurrent$contactRates[modelAgeGroupCatIndex]
	rhoAge <- pCurrent$rho[modelAgeGroupCatIndex] 
	
	
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
	
	psi <- pCurrent$R0*hostSurvivalTotal*pCurrent$LDecayRate/(pCurrent$lambda*pCurrent$z*summation)
	
	return(psi*2) ## This two because it's infection with both male and female worms and only half of these are going to be female. 
}

## new psi calculation. 
## Results are returned as a list because we want deltaT and B for IC calculation. 
getPsi <- function(pCurrent)
{
	## higher resolution. 
	deltaT <- 0.1
	
	## inteval-centered ages for the age intervals. 
	## 0 to maxHostAge, mid points. 
	modelAges <- seq(0,pCurrent$maxHostAge-deltaT,by=deltaT) + 0.5*deltaT ## each age group in the actual model is annual characterised by mid-value (0.5 is [0,1), etc.)
	
	## hostMu for the new age intervals. 
	hostMuGroupIndex <- cut(modelAges,breaks=pCurrent$muBreaks,labels=1:length(pCurrent$hostMuData)) 
	hostMu <- pCurrent$hostMuData[hostMuGroupIndex] 
	
	meanDeaths <- hostMu*deltaT
	hostSurvivalCurve <- exp(-cumsum(meanDeaths))
	
	MeanLifespan <- sum(hostSurvivalCurve[1:length(modelAges)])*deltaT  ## This should be the integral under the survival curve UP TO THE TOP AGE LIMIT. 
	
	## calculate the cumulative sum of host and worm death rates from which to calculate worm survival. 
	intMeanWormDeathEvents <- cumsum(hostMu+pCurrent$sigma)*deltaT
	
	## need rho and beta at this age resolution as well. 
	ageIndices <- 1:(length(pCurrent$contactAgeGroupBreaks)-1) 
	modelAgeGroupCatIndex <- cut(modelAges,breaks=pCurrent$contactAgeGroupBreaks,labels=ageIndices) 
	betaAge <- pCurrent$contactRates[modelAgeGroupCatIndex]
	rhoAge <- pCurrent$rho[modelAgeGroupCatIndex] 
	
	
	wSurvival <- exp(-pCurrent$sigma*modelAges) 
	
	## calculate the infectiousness bit first. 
	nMax <- length(hostMu)
	B <- rep(0,nMax)
	
	for(i in 1:nMax)
	{
		B[i] <- sum(betaAge[1:i]*wSurvival[i:1])*deltaT
	}
	
	summation <- sum(rhoAge*hostSurvivalCurve*B)*deltaT
	
	psi <- pCurrent$R0*MeanLifespan*pCurrent$LDecayRate/(pCurrent$lambda*pCurrent$z*summation)
		
	return(psi)  
}


# Draw a lifespan from the population survival curve. 
getLifeSpans <- function(nSpans,pars)
{
	spans <- rep(0,nSpans)
	for(i in 1:nSpans)
	{
		spans[i] <- which(runif(1)*max(pars$hostAgeCumulDistr)<pars$hostAgeCumulDistr)[1]
	}
	return(pars$muAges[spans])
}



## takes the parameter structure. 
## returns a structure containing the equilibrium worm burden with age and the reservoir value. 
## Also the breakpoint reservoir value and other things. 
## N.B. - psi should already be calculated. 
getEquilibrium <- function(pCurrent)
{
  ## higher resolution. 
  deltaT <- 0.1
  
  ## inteval-centered ages for the age intervals. 
  modelAges <- seq(0,pCurrent$maxHostAge-deltaT,by=deltaT) + 0.5*deltaT ## each age group in the actual model is annual characterised by mid-value (0.5 is [0,1), etc.)
  
  ## hostMu for the new age intervals. 
  hostMuGroupIndex <- cut(modelAges,breaks=pCurrent$muBreaks,labels=1:length(pCurrent$hostMuData)) 
  hostMu <- pCurrent$hostMuData[hostMuGroupIndex] 
  
  meanDeaths <- hostMu*deltaT
  hostSurvivalCurve <- exp(-cumsum(meanDeaths))
  
  MeanLifespan <- sum(hostSurvivalCurve[1:length(modelAges)])*deltaT  ## This should be the integral under the survival curve UP TO THE TOP AGE LIMIT. 
  
  ## need rho and beta at this age resolution as well. 
  ageIndices <- 1:(length(pCurrent$contactAgeBreaks)-1) 
  modelAgeGroupCatIndex <- cut(modelAges,breaks=pCurrent$contactAgeBreaks,labels=ageIndices) 
  betaAge <- pCurrent$contactRates[modelAgeGroupCatIndex]
  rhoAge <- pCurrent$rho[modelAgeGroupCatIndex] 
  
  wSurvival <- exp(-pCurrent$sigma*modelAges) 
  
  nMax <- length(hostMu)
  Q <- rep(0,nMax)  ## This variable times L is the equilibrium worm burden. 
  
  for(i in 1:nMax)
  {
    Q[i] <- sum(betaAge[1:i]*wSurvival[i:1])*deltaT
  }
  
  ## converts L values into mean force of infection. 
  FOIMultiplier <- sum(betaAge*hostSurvivalCurve)*deltaT/MeanLifespan
  
  ## upper bound on L. 
  SRhoT <- sum(hostSurvivalCurve*rhoAge)*deltaT
  R_power <- 1/(pCurrent$k+1)
  L_hat <- pCurrent$z*pCurrent$lambda*pCurrent$psi*SRhoT*pCurrent$k*(pCurrent$R0^R_power - 1)/(pCurrent$R0*MeanLifespan*pCurrent$LDecayRate*(1-pCurrent$z)) 
  
  ## now check the value of the K function across a series of L values.
  ## find point near breakpoint. L_minus is the value that gives an age-averaged worm burden of 1. Negative growth should exist somewhere below this. 
  L_minus <- MeanLifespan/sum(Q*hostSurvivalCurve*deltaT)
  test_L <- c(seq(0,L_minus,length.out=10),seq(L_minus,L_hat,length.out=20))
  
  K_valueFunc <- function(currentL,pCurrent)
  {
    answer <- pCurrent$psi*sum(pCurrent$reproFunc(currentL*Q,pCurrent)*rhoAge*hostSurvivalCurve*deltaT)/(MeanLifespan*pCurrent$LDecayRate) - currentL
    return(answer)
  }
  
  K_values <- sapply(test_L, FUN=K_valueFunc,pCurrent=pCurrent)
  
  ## now find the maximum of K_values and use bisection to find critical Ls.
  iMax <- which.max(K_values)
  
  mid_L <- test_L[iMax]
  
  ## is mid_l < 0? 
  ## 20/5/16: why not just check the K value at this point? This also guarrantees that the uniroot call below has interval ends of opposite sign.  
  if(K_values[iMax] < 0) 
  {
    solutions <- list(stableProfile=0*Q,ageValues=modelAges,L_stable=0,L_breakpoint=NA,K_values=K_values,L_values=test_L,FOIMultiplier=FOIMultiplier)
    return(solutions)
  }
  
  ## find the top L...
  L_stable <- uniroot(K_valueFunc,interval=c(mid_L,4*L_hat), extendInt="yes", pCurrent=pCurrent)$root
  
  ## find the unstable L... Start at 1 in from the zero at the bottom. 
  #L_break <- NA  ## sometimes the first point is not within the negative range. Need to use worms 0 to 2 as defined on the L scale.  
  L_break <- test_L[2]/50 
  if(K_valueFunc(L_break,pCurrent)<0)  ## if it is less than zero at this point, find the zero. 
  {
    L_break <- uniroot(K_valueFunc,interval=c(L_break,mid_L), extendInt="yes", pCurrent=pCurrent)$root
  }
  
  
  stableProfile <- L_stable*Q 
  hatProfile <- L_hat*Q
  
  solutions <- list(stableProfile=stableProfile, hatProfile=hatProfile, ageValues=modelAges, hostSurvival=hostSurvivalCurve, L_breakpoint=L_break, K_values=K_values, L_values=test_L, FOIMultiplier=FOIMultiplier, L_stable=L_stable, L_hat=L_hat)
  return(solutions)
}





