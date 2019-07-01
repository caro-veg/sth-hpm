######################################################################################################
#
# Results format is a list containing the parameters ($params, last list element) 
# and results from individual realisations (list elements 1,2,3...)
#
# The results from each individual realization is a list. 
# Its elements are:
# lists for each time point specified by the outTimings parameter. 
# $repNo - the realization number (don't have a use for this at the moment). 
#
# The list for each time point in each realization contains: 
# worms: a data.frame of total and female worms at that time for each host.
# hosts: a data.frame with birth dates and death dates for each host
# freeLiving: reservoir
# time: time at which the data stored in this list element was outputed
# adherenceFactors: probabilities with which each host takes their treatment
#
#######################################################################################################



# Remove objects
rm(list=ls())

closeAllConnections()
runtime = proc.time()


############################################################################
# WHEN RUNNING SIMULATIONS ON CLUSTER USE THIS BLOCK TO READ IN 
# PARAMETERS FROM COMMAND LINE/BATCH FILE
############################################################################

## command line arguments. 
args <- commandArgs(trailingOnly = TRUE)

outFileStub <- args[1]				# stub for output file name
dataInputPath <- args[2]			# path to parameter files
outPath <- args[3]				# path to where data should be outputted
paramFileName <- args[4]			# name of file with simulation parameters
currentDemogName <- args[5]			# name of demographies used in simulation
							# read from separate file

numberPeopleToMigrate <- as.numeric(args[6])	# number of people to migrate between villages

seed <- as.numeric(args[7])				# seed for repeated runs


############################################################################
# SET WORKING DIRECTORY AND PARAMETER FILE
############################################################################

setwd(dataInputPath)						# function files etc should be in input directory
parameterPath <- paste0(dataInputPath, paramFileName)	# path plus filename of 


############################################################################
# SET UP INITIAL CONDITIONS FOR SIMULATION
############################################################################

# Source functions for stochastic code
source(paste0(dataInputPath, "ParasiteFunctions.R"))
source(paste0(dataInputPath, "helsim_FUNC.R"))


# Set up the output and recording
#logFile <- file(paste0(outPath,outFileStub,"_log.txt"), "w")
logFile <- file(paste0(outPath, outFileStub, "_", numberPeopleToMigrate, "_log.txt"), "w")
cat(file=logFile,"Log file is open!\n")

# Read in parameters
testExist <- file(parameterPath,"r")
if(!isOpen(testExist)){
	cat(file=logFile,"Parameter file not found: ",parameterPath,"\n")
	stop("Parameter file not found.")
} else {
	cat(file=logFile,"Parameter file found: ",parameterPath,"\n")
}
close(testExist)
params <- readParams(fileName=parameterPath, demogName=currentDemogName)


# Configure parameters
params <- configure(params)
params$psi <- getPsi(params)  	 


# Call function to return equilibrium reservoir level for IC. 
equi <- getEquilibrium(params)
params$equiData <- equi


params$numberPeopleToMigrate <- numberPeopleToMigrate
params$seed <- seed

###################################################################


# Specify a single realisation/rep as a function
doRealization <- function(params)
{
	# Setup the simulation
	simData <- setupSD(params)	# Set up simulation data for village all worm counts and free living worms are set to 0 below
	simData_static <- simData			# Copy original village as future source of infection
			
	simData$freeLiving <- 0
	simData$worms$total <- rep(0, params$N)
	simData$worms$female <- rep(0, params$N)

	

	time <- 0							# Set initial time to 0
	FLlast <- time						# Time at which to update simulation data using the doFreeLive function
	outTimes <- params$outTimings				# Times at which data should be recorded
	
	emigrationTimes <- params$emigrationTimes		# Times at which people leave their home villages
	immigrationTimes <- params$immigrationTimes	# Times at which people return to their home villages
	
	nextOutIndex <- which.min(outTimes)								# Determine time when data should be recorded next
	nextOutTime <- outTimes[nextOutIndex]
	ageingInt <- 1/52     										# Check ages every week - should this be in a parameter file?
	nextAgeTime <- ageingInt									# Time at which individuals' age is advanced next
	nextEmigrationIndex <- which.min(emigrationTimes)
	nextEmigrationTime <- emigrationTimes[nextEmigrationIndex]			# Time at which individuals emigrate next
	nextImmigrationIndex <- which.min(immigrationTimes)
	nextImmigrationTime <- immigrationTimes[nextImmigrationIndex]		# Time at which individuals immigrate next (should always be later than next emigration in seasonal scenario)
	
	#print(paste("nextEmigrationTime:", nextEmigrationTime))
	#print(paste("nextImmigrationTime:", nextImmigrationTime))


	nextStep <- min(nextOutTime, time+params$maxStep, nextAgeTime, nextEmigrationTime, nextImmigrationTime)	# Determine which event will take place next
	breakTime <- NA 	# Set time at which breakpoint has been crossed to NA because at beginning of simulation this has not happened yet
				# will be set to time when breakpoint has been crossed in simulation
	
	results <- list() # initialise empty list to store results 	
	outCount <- 1	# keep count of outputs

	migrantIndices <- c()	# vector to store people that will migrate
	
	migrantAgeGroupIndices <- which((time - simData$demography$birthDate) >= params$migrantAgeGroup[1] & (time - simData$demography$birthDate) < params$migrantAgeGroup[2])			
	if(length(migrantAgeGroupIndices) > params$numberPeopleToMigrate)
	{
		migrantIndices <- sample(migrantAgeGroupIndices, params$numberPeopleToMigrate)
		simData$worms$away[migrantIndices] <- TRUE
	}else
	{
		migrantIndices <- migrantAgeGroupIndices
		simData$worms$away[migrantIndices] <- TRUE
	}


	
	# Run stochastic algorithm	
	while((time<params$maxTime))    
	{
		
		rates <- calcRates(params, simData, time)
		sumRates <- sum(rates)

		# Determine if worms have gone above breakpoint in village 
		if((simData$freeLiving>params$equiData$L_breakpoint) & is.na(breakTime))  
		{
			breakTime <- time
			cat(file=logFile, "Total worms: ", simData$worms$total, "\n") 
			cat(file=logFile, "breakTime: ", breakTime, "\n")
		}
    
		if(sumRates > 0.001) 
		{
			tstep <- rexp(1, sumRates)	# draw time interval from exponential distribution
		}else
		{
			tstep <- params$maxTime
		}

		print(paste("sumRates:", sumRates))
		print(paste("time:", time))
		print(paste("tstep:", tstep))
		print(paste("nextStep:", nextStep))

				
		if((time+tstep)<nextStep) 
		{
			time <- time + tstep			
			simData <- doEvent(rates, params, simData, time)

		}else 
		{ 
			time <- nextStep

		  	## Predetermined event block.
			simData <- doFreeLiveMigration(params, simData, nextStep-FLlast)
			FLlast <- nextStep			 
			timeBarrier <- nextStep + 0.001

			# Ageing and death 
			if(timeBarrier>nextAgeTime)
			{ 
		  		simData <- doDeathSimple(params, simData, time) 
				nextAgeTime <- nextAgeTime + ageingInt
			}

			
	
			# EMIGRATION: set away flag to true for given number of migrants
			if(timeBarrier>nextEmigrationTime)
			{	
				emigrationTimes[nextEmigrationIndex] <- params$maxTime + 10
				nextEmigrationIndex <- which.min(emigrationTimes)
				nextEmigrationTime <- emigrationTimes[nextEmigrationIndex]
				print(paste("nextEmigrationTime:", nextEmigrationTime))

				print(migrantIndices)
			}
			
			
			# IMMIGRATION: At a  specified time people return into their home village 
			# before they return they get infected 
			# sample worm burden from equilibrium age structure and worm burden from 
			if(timeBarrier>nextImmigrationTime)
			{
				infectantAgeGroupIndices <- which((time - simData_static$demography$birthDate) >= params$migrantAgeGroup[1] & (time - simData_static$demography$birthDate) < params$migrantAgeGroup[2])
				if(length(infectantAgeGroupIndices) > length(migrantIndices))
				{
					infectantIndices <- sample(infectantAgeGroupIndices, length(migrantIndices))
					simData$worms$away[migrantIndices] <- TRUE
				}else
				{
					infectantIndices <- sample(1:params$N, length(migrantIndices))
					simData$worms$away[migrantIndices] <- TRUE
				}				

				simData$worms[migrantIndices, "total"] <- simData_static$worms[infectantIndices, "total"]
				simData$worms[migrantIndices, "female"] <- simData_static$worms[infectantIndices, "female"]
				simData$worms[migrantIndices, "away"] <- FALSE

				# Chemotherapy only for returning migrants - only if coverage > 0
				if(params$coverage > 0) 
				{
					simData <- doChemoSimpleImmigrantsOnly(params, simData, migrantIndices, time)
				}
				
				immigrationTimes[nextImmigrationIndex] <- params$maxTime + 10
				nextImmigrationIndex <- which.min(immigrationTimes)
				nextImmigrationTime <- immigrationTimes[nextImmigrationIndex]
				print(paste("nextImmigrationTime:", nextImmigrationTime))
				print(migrantIndices)
				print(sum(simData$worms[migrantIndices, "total"]))
			}	

			
	
			# Record
			if(timeBarrier>nextOutTime )  
			{ 	
				currentRes <- list()
				currentRes$worms <- simData$worms
				currentRes$hosts <- simData$demography
				currentRes$freeLiving <- simData$freeLiving
				currentRes$time <- time

				results[[outCount]] <- currentRes
				outCount <- outCount+1
				
				outTimes[nextOutIndex] <- params$maxTime+10
				nextOutIndex <- which.min(outTimes)
				nextOutTime <- outTimes[nextOutIndex]
			}		
			
			nextStep <- min(nextOutTime, time+params$maxStep, nextAgeTime)
			
		} ## end of predetermined event block.
	} ## end of while loop. 
	
	
	finalResults <- list(results=results, breakTime=breakTime)
	
	return(finalResults)
}


cat(file=logFile, "Setup complete after ", (proc.time() - runtime)[["elapsed"]], " secs.\n")


## Set up and run multiple instances of doRealization() in parallel

cat(file=logFile, "Setting up parallel stuff.\n")

require("foreach")
require("doParallel", quietly=TRUE)
require("doRNG", quietly=TRUE)

cl <- makeCluster(params$nNodes)
registerDoParallel(cl)
registerDoRNG(seed=params$seed)

cat(file=logFile,"Starting parallel stuff.\n")

foreachResults <- foreach(i=1:params$numReps) %dorng% doRealization(params)

stopImplicitCluster() # On advice from Wes. 

foreachResults$params <- params


#cat(file=logFile,"Total elapsed time = ",(proc.time() - runtime)[["elapsed"]],"\n")


resultsFile <- paste0(outPath, outFileStub, "results_as_", params$numberPeopleToMigrate, "ppl_", "R0_", params$R0, "_k", params$k, "_seed", params$seed, ".RData")


save(file=resultsFile,results=foreachResults)


cat(file=logFile,"Log file closing...\n")
close(logFile)



