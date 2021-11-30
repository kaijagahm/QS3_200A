
#### Vellend community simulation code ####

## Original document by:
## Nathan Kraft, UCLA
## nkraft@ucla.edu

## Updated Nov. 2021 by:
## Colin Kremer, UCLA
## kremer@ucla.edu

####################################

## Derived from code by Mark Vellend and Andrew McDonald, from Vellend's Theory of Ecological Communities (Princeton, 2016)
## basic simulation of a two species explicit-metacommunity neutral model with migration but no speciation, 
## adding the potential to introduce fitness differences and frequency dependence that varies among local communities. 


####################################
## main simulation as a function: #
##################################

## first function simulates 2 species in metacommunity, tracking frequency of species 1 over time in output:

simulate_metacommunity <- function(years = 50, patches = 10, J = 100, m = 0,
                                   fitness_ratio_ave = rep(1, patches), 
                                   frequency_dependence = rep(0, patches)){
	
	###############################
	## Description of Arguments: #
	#############################
	
	## years - here defined arbitrarily as (J * patches) deaths = 1 year
	## this is helpful if you want to vary J and compare results on same timescale
	
	## patches - number of local communities to track
	
	## J - the number of individuals per patch
	
	## m - migration rate, probability that a local birth
	## comes from the metacommunity instead of reproduction of a local individual
	## The value of m needs to be between 0 and 1. 
	
	
	## fitness_ratio_ave - specifies the fitness difference between species 1 and 2 in each patch.
	## Needs to be vector of length 'patches'; the default value is 1 (neutral)
	## GUIDANCE ON VALUES:
	## 1 = no fitness difference between species 1 and 2 (neutral)
	## If you want symmetric fitness differences between the two species in different patches,
	## the your values need to be INVERSES of each other, such as c(1.2, 1/1.2), NOT c(1.2, 0.8)
	
	
	## frequency_dependence - specifies the slope of frequency dependence in each patch.
	## Needs to be vector of length 'patches'; the default value is 0 (no frequency dependence)
	## 0 implies no frequency dependence. 
	## > 0 implies positive frequency dependence (works in favor of abundant species)
	## < 0 implies negative frequency dependence (works against abundant species)
	## Note that fitness_ratio_average above determines the intercept of the 
	## relationship for frequency dependence 
	
	### a few fail-safe tests to stop the simulation if key arguments are outside of allowed ranges:
	if(m<0 | m>1){return(print("simulation canceled: migration m needs to be between 0 and 1"))}
	if(length(fitness_ratio_ave)!=patches){return(print("simulation canceled: fitness_ratio_ave needs to be of length 'patches' "))}
	if(length(frequency_dependence)!=patches){return(print("simulation canceled: fitness_ratio_ave needs to be of length 'patches' "))}
	
	########################
	### start simulation ##
	######################
	
  ## Here the metacommunity is modeled explicitly, so Jm (total metacommunity size) is:
  Jm <- J * patches
  
	## Initialize output matrix, where rows are years and columns are local communities (patches). 
	## We'll be working with two species in a zero sum model, so we only need to track one species
	frequency_sp1 <- matrix(nrow = years, ncol = patches) 
	
	## define initial population size of species 1. 
	## This is set by default to be equal in abundance to species 2:
	initial_abundance_sp1 <- 0.5*J 
	
	## our live, working snapshot of the metacommunity at one point in time;
	## rows are individuals and each column is one local community (or patch). 
	metacommunity_state <- matrix(nrow=J, ncol=patches)
	
	## set half the individuals in all patches to species 1:
	metacommunity_state[1:initial_abundance_sp1,] <- 1
	## and the other half to species 2:
	metacommunity_state[(initial_abundance_sp1+1):J,] <- 2 
	
	
	## record frequency of species 1 for year 1 - this is same for all patches initially.
	frequency_sp1[1,]<-initial_abundance_sp1/J
	
	## for clarity here I will nest 'for' loops here for years and deaths, though this makes it slower to run:
	
	## loop over years:
	for(i in 2:years){
		## assuming our initial data represents year 1, we'll start at year 2:
		
	  ## loop over all individuals
		for(j in 1:Jm){
			## each year has J * patches (Jm) deaths in it- we'll loop over that here:
			
			## for each death, choose one patch as location- all have same chance:
			patch<-sample(1:patches, 1)
			
			## calculate probability the dead individual is replaced by species 1:
			
			if(runif(1)<m){
				## randomly choose a number between 0 and 1,
				## if it is less than m, dead individual is replaced by
				## dispersal from metacommunity.
				## chance of replacement from the metacommunity being species 1 is just 
				## the relative abundance of species 1 in species pool
				Prob_repro_Sp1<-sum(metacommunity_state==1)/Jm
			} else {
				## if the random number was not less than m, the dead individual will be
			  ## replaced from the local community (patch).
				## If present, fitness differences and frequency dependence
				## play out in this local replacement process
				
				# calculate frequencies of species 1 and 2 in selected patch:
				local_freq_1 <- sum(metacommunity_state[,patch]==1)/J
				local_freq_2 <- 1 - local_freq_1  # if we know the frequency of species 1, we can calculate that of species 2
				
				## Then, calculate the fitness ratio between species 1 and 2 in this randomly selected patch.
				## This value will depend on the frequency dependence for this particular patch 
				## (taken from the vector provided to the function)
				## It also depends on the fitness ratio average for the particular patch (also taken from the provided vector)
				
				fitness_ratio <- fitness_ratio_ave[patch]*exp(frequency_dependence[patch]*(local_freq_1-0.5))
				#fitness_ratio <- exp(frequency_dependence[patch]*(local_freq_1-0.5)+log(fitness_ratio_ave[patch]))
				
				##  use this ratio and the observed frequencies of species to calculate chance that species 1 reproduces:
				Prob_repro_Sp1 <- fitness_ratio*local_freq_1 / ( fitness_ratio*local_freq_1 + local_freq_2)
				
				## if this ratio is > 1, the chance that species 1 reproduces is elevated relative to its abundance.
				## if this ratio is < 1, the chance species 1 reproduces is depressed relative to its abundance.
				
			}
			
			## now randomly choose the individual inside the selected patch that will die
			unlucky.individual <- ceiling(J*runif(1))
			
			## Replace this unlucky individual by a new individual of species 1 with the probability calculated above,
			## using the sample() command
			metacommunity_state[unlucky.individual, patch] <- sample(c(1,2), 1, prob=c(Prob_repro_Sp1, 1 - Prob_repro_Sp1))
			
		}
		
		## at the end of each "year" (= after Jm deaths have occurred), record the state of the metacommunity:
		frequency_sp1[i,] <-colSums(metacommunity_state==1)/J ## counts frequency of sp 1 in each patch
		
	}
	
	## after the years are done, return the final data set documenting the frequency of sp 1 over time in each patch:
	return(frequency_sp1)
}



#########################################################################
## function for simple line plots of output from above simulation: ###
####################################################################


trace_metacommunity=function(df, title="your title here",colorPatches=FALSE){
	## expects output from simulate_metacommunity()
	## you can add any title you like to track parameters, etc
	
  if(colorPatches){
    colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
    cols <- colfunc(ncol(df))
  }else{
    cols <- rep('black',ncol(df))
  }
  
  plot(1:nrow(df), df[,1], type='l', xlab="years", ylab="frequency of species 1", 
       ylim=c(0,1), main=title, col=cols[1])
  for(i in 2:(ncol(df))){
    lines(1:nrow(df), df[,i], type='l',col=cols[i])
  }
}


###################################
## simple examples ####
###################################

## run functions above first to save them into memory.

# PROBLEM 2
# Using a relatively small local community size (J=100), number of patches (10), and no selection, explore how variation in the strength of migration interacts with ecological drift to shape variation in species composition across communities (patches). In multi-species communities we would typically use a measure of beta diversity, but as we only have two species here, please use the variance in the relative abundance of species 1 at 50 years as your measure of heterogeneity in composition across communities (patches). Produce a plot (or plots) that demonstrates your findings, and describe the overall interaction that you see between these two processes. Include the R code you used to produce the figure.

m_values <- round(seq(from = 0, to = 1, length.out = 10), 2)

variances <- rep(NA, length(m_values))

pdf("question2.pdf", width=8, height=9)
layout(matrix(1:10, nrow=5, byrow=T))
for(i in 1:length(m_values)){
  title <- paste0("neutral, J = 100, m=", m_values[i])
  output <- simulate_metacommunity(years = 50, patches = 10, J = 100, m = m_values[i], fitness_ratio_ave = rep(1, 10), frequency_dependence = rep(0,10))
  trace_metacommunity(output, title, colorPatches = F)
  
  # calculate variance in frequency of species 1 at time = 50 years
  heterogeneity <- var(output[50,]) # variance of the 50th row
  variances[i] <- heterogeneity
  
  print(paste0("finished iteration ", i))
}
dev.off()

# We can see from these graphs that as m increases toward 1 (i.e. there is more connectivity/flow of individuals from other patches), the frequency of species 1 stays more consistent. Therefore, connectivity between patches in a metapopulation under these conditions (i.e. no selection) facilitates coexistence between the two species. By contrast, when m is very small, there's a lot of variability. Sometimes species 1 completely takes over, and sometimes it goes totally extinct (i.e. species 2 takes over).

# Let's plot the variances at year 50 for different values of m
plot(x = m_values, y = variances, type = "l", xlab = "migration rate", ylab = "Variance in sp1 freq at year 50")

# Interestingly, this relationship is not linear. When m is 0, the variance in the frequency of species 1 at 50 years is very high, but with even a small amount of migration, the variance drops way, way down. Overall, it seems that migration counters drift.

# PROBLEM 3
# Explore the interaction of migration and spatially variable but constant selection, such as you might expect with habitat filtering across a landscape. I suggest starting with a relatively small community size for feasibility (J=100, it's OK if there is some drift happening in the background). Specifically, explore what happens when your communities are equally divided between two habitat types, one type that favors species one and another type that favors species two. The fitness advantages conferred by the two habitats should be symmetric (the fitness_ratio_ave values in each habitat should be the inverse of one another), and there should be no frequency dependence. How does increasing migration influence variation in community composition across the landscape? How does diversity within local communities vary with increasing migration? Summarize your understanding of how spatially variable, constant selection and migration interact. If you see evidence of any ecological theories that you know from elsewhere represented in these results, please name them. Include any code needed to replicate your findings.

# Each patch is either habitat type 1 (favors species 1) or habitat type 2 (favors species 2).

variances <- rep(NA, length(m_values))

pdf("question3.pdf", width=8, height=9)
layout(matrix(1:10, nrow=5, byrow=T))
for(i in 1:length(m_values)){
  title <- paste0("const. sel., 2 habs, fitRatios = 1.2 & 1/1.2, m=", m_values[i])
  output <- simulate_metacommunity(years = 50, patches = 10, J = 100, m = m_values[i], fitness_ratio_ave=c(rep(1.2,5), rep(1/1.2, 5)), frequency_dependence = rep(0,10))
  trace_metacommunity(output, title, colorPatches = F)
  
  # calculate variance in frequency of species 1 at time = 50 years
  heterogeneity <- var(output[50,]) # variance of the 50th row
  variances[i] <- heterogeneity
  
  print(paste0("finished iteration ", i))
}
dev.off()

# Once again, greater migration rates facilitate coexistence. If there's no migration, it is absolute that over time, species 1 will completely dominate the patches where it has an advantage, while species 2 will do the same in the other patches (in this case, because there are five of each habitat type, each species will take over in five of the patches).

# Let's plot the variances at year 50 for different values of m
plot(x = m_values, y = variances, type = "l", xlab = "migration rate", ylab = "Variance in sp1 freq at year 50")

# In contrast to the previous graph, here we don't see a large reduction in the 50 year variance until m is slightly higher. This indicates that when selection is happening, you need more migration to compensate for those effects than you would when the habitats are all equally suitable to both species. 

# This reminds us of the idea of source-sink dynamics--even if a patch isn't suitable to a given species, that species might sometimes persist by continually migrating in from other, more suitable habitats. 

# QUESTION 4
# Investigator's choice! Pick a third combination of two or more processes and explore how they interact. Begin by stating your question and describing the simulations you will use to address this. Summarize your findings graphically and interpret them. Include any code needed to replicate your findings.

# How does migration affect the populations when there is negative frequency dependence to varying extents (i.e. when a species is abundant, it is at a disadvantage)?

m_values <- round(seq(from = 0, to = 1, length.out = 4), 2)
freq_dependence_values <- c(0, -0.5)

pdf("question4.pdf", width=8, height=9)
layout(matrix(1:8, nrow=4, byrow=F))

for(i in 1:length(m_values)){
  for(j in 1:length(freq_dependence_values)){
    title <- paste0("m = ", m_values[i], "; freq dep = ", freq_dependence_values[j])
    output <- simulate_metacommunity(years = 50, patches = 10, J = 100, m = m_values[i],
                                     fitness_ratio_ave = rep(1, 10), frequency_dependence = rep(freq_dependence_values[j], 10))
    trace_metacommunity(output, title, colorPatches = F)
    
    print(paste0("completed iteration ", j, " for m_value ", i))
  }
}

dev.off()

# We can see from looking at these plots that both migration rate and negative frequency dependence promote coexistence (i.e. increase stability of the frequency of species 1). It looks to us like migration rate has a larger effect on stability than negative frequency dependence does, but comparing within pairs of situations when m is held constant but frequency dependence is changed, we do see a clear reduction in variability when there is  greater negative frequency dependence.

# This reminds us of Lotka-Volterra dynamics, where coexistence is promoted by each species having density-dependence, i.e. being more limited by itself than by the other species.
