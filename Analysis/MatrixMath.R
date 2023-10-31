#' @description 
#' Initiates a random fun of spatially-explicit trait matching model given an initial set number plant and animal species.
#' Trait values are randomly assigned, and new species are included through time. Outputs results of simulation in a list object.
#' 
#' @param disptype character; describes relationship between population growth rate and dispersal strength. Options are one of c("neutralComp", "positiveComp", "negativeComp")
#' When selecting "positiveComp", animal populations with disperse at rate lambda/2 \propto disprobmax (with a maximum of disprobmax). 
#' "negativeComp" reverses this relationship, such that animal pops disperse at rate \propto 1-lambda/2 disprobmax (same limits as above)\
#' "neutralComp" removes the correlation, setting animal dispersal equal to disprobmax/2 for all species.
#' @param nsites number; Number of sites in metacommunity
#' @param dexpsim Distance matrix; describing distance between sites
#' @param disprobmax
#' @param n_plants integer; initial number of plant species in community
#' @param n_animals integer; initial number of animal species in community
#' @param r number (0-infinity); maximum per capita growth rate for plant species
#' @param mup number (0-1); plant death rate per timestep
#' @param mua number (0-1); animal death rate per timestep
#' @param lambda number (0-infinity); maximum per capita growth rate for animal species
#' @param K integer; carrying capacity for plant species (same for all species)
#' @param e_thresh integer; extinction threshold; species with total populations across sites below this threshold are considered extinct 
#' @param invProb number (0-1); probability on new species invasion per timestep. Average time between invasions is equal to 1/invProb.
#' Once invasion occurs, we flip a coin to determine whether the new species is a plant or animal.
#' @param invade_size integer; size of invading population. All invaders of a given species colonize a randomly selected site together
#' @param num_timeSteps integer; number of timesteps to run the simulation for (default value is 500)
#' 
#' 
#' @return outputlist list object; with the following levels;
#' 
#' @return outputlist$p_pops_output matrix of integers; first column is site number, second column is time, columns 3:ncol are population sizes of plant_species
#' @return outputlist$a_pops_output matrix of integers; analogous form as above
#' @return outputlist$a_traitM vector of numbers; mean values for each animal species' trait distributions
#' @return outputlist$a_traitV vector of numbers; variance parameter for each animal species' trait distributions
#' @return outputlist$p_traitsM vector of numbers; analogous form as above
#' @return outputlist$p_traitV vector of numbers; analogous form as above
#' 

runDispersalSim <- function(X, disptype, nsites, dexpsim, disprobmax, n_plants, n_animals, r, mup, mua, o, lambda, K, e_thresh, invProb, invade_size, num_timeSteps=500, seedstart=0){
  factor_sites <- as.factor(1:nsites)
  #Setting up initial plant pops
  p_pops <- matrix(data=round(runif(nsites*n_plants, 1, 200)), nrow=nsites, ncol=n_plants) #now, our pops are a matrix instead of a vector in order to track across multiple sites
  p_traitM <- runif(n_plants, 0, 1)
  p_traitV <- runif(n_plants, 0, .25)
  #p_traitM <- c(0.25, 0.99, 0.5)
  #p_traitV <- c(0.1, 0.1, 0.1) #All set to be the same in Becker paper
  
  #Setting up initial animal pops
  a_pops <- matrix(data=round(runif(nsites*n_animals, 1, 200)), nrow=nsites, ncol=n_animals)
  a_traitM <- runif(n_animals, 0.05, 1)
  a_traitV <- runif(n_animals, 0.05, .25)
  
  #Add our function for finding overlap
  min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
    f1 <- dnorm(x, mean=mu1, sd=sd1)
    f2 <- dnorm(x, mean=mu2, sd=sd2)
    pmin(f1, f2)
  }
  
  #pops <- NULL
  p_pops_output <- cbind(c(1:nsites), rep(0,nsites),p_pops) #Make output dataframe; first column is site, second is timestep; then populations
  a_pops_output <- cbind(c(1:nsites), rep(0,nsites),a_pops) #Make output dataframe; first column is site, second is timestep; then populations

  prior_richness <- 0 #setting as the starting point
  
  int_fail <- FALSE #Dummy variable for integral try statement
  lambda_fail <- FALSE #Dummy variable for rpois() projection
  
  P_TraitVTot <- NULL #Create a bunch of empty objects so we can rbind weighted mean trait variance to later
  A_TraitVTot <- NULL
  P_WTraitMean_output <- NULL
  A_WTraitMean_output <-  NULL
  
  propOutputtemp <- list() #Set up our realized specialism output object
  pollBenefitsOutput <- list() #same for plants
  
  for(t in 1:num_timeSteps){
    set.seed(t+seedstart)
    propOutputtemp[[t]] <- list() #Make blank list entry for this timestep that we can fill in later with our plant benefit
    pollBenefitsOutput[[t]] <- list() #Same thing above for the plants
    richness <- ncol(a_pops)+ncol(p_pops)
    if(richness > prior_richness){ #Check if we've updated the number of spp in the community
      #if so, recalculate niche overlap and competition matrices
      #print("recalculating")
        alpha <- matrix(nrow=n_animals, ncol=n_plants) #Create an empty matrix for plant-pollinator partnerships
        for(ia in 1:n_animals){
          if(int_fail==TRUE){ #If we've come across an integration failure, break the loop
            break()
          }
          for(ip in 1:n_plants){
            temp <- try(integrate(min.f1f2, -Inf, Inf, mu1=p_traitM[ip], mu2=a_traitM[ia], sd1=p_traitV[ip], sd2=a_traitV[ia])$value) #Integrate our minimum equation over all numbers; output is the total area overlapping both curves, with max value of perfectly overlapping curves as 1
            if(class(temp)=='try-error'){
              int_fail <- TRUE
              break()
            }
            alpha[ia, ip] <- temp
            #print(paste("ia=", ia, "; ip=", ip, sep=""))
          }
        }
        alpha[alpha<0.05] <- 0
        
        u <- matrix(data=NA, nrow=n_plants, ncol=n_plants)
        
        for(p1 in 1:n_plants){
          if(int_fail==TRUE){ #If we've come across an integration failure, break the loop
            break()
          }
          for(p2 in 1:n_plants){
            temp <- try(integrate(min.f1f2, -Inf, Inf, mu1=p_traitM[p1], mu2=p_traitM[p2], sd1=p_traitV[p1], sd2=p_traitV[p2])$value) #Integrate our minimum equation over all numbers; output is the total area overlapping both curves, with max value of perfectly overlapping curves as 1
            if(class(temp)=='try-error'){
              int_fail <- TRUE
              break()
            }
            u[p1, p2] <- temp
          }
        }
    
    o <- matrix(data=NA, nrow=n_animals, ncol=n_animals) #Now, let's look at animal competition
    for(ia in 1:n_animals){
      if(int_fail==TRUE){ #If we've come across an integration failure, break the loop
        break()
      }
      for(a2 in 1:n_animals){
        temp <- try(integrate(min.f1f2, -Inf, Inf, mu1=a_traitM[a2], mu2=a_traitM[ia], sd1=a_traitV[a2], sd2=a_traitV[ia])$value) #Integrate our minimum equation over all numbers; output is the total area overlapping both curves, with max value of perfectly overlapping curves as 1
        if(class(temp)=='try-error'){
          int_fail <- TRUE
          break()
        }
        o[ia, a2]<- temp
        #print(paste("ia=", ia, "; ip=", ip, sep=""))
      }
    }
  }
    
    prior_richness <- richness
    
    if(int_fail==TRUE){ #Same check for int_fail as above, but this time we'd be breaking our overall loop
      outputlist <- list(plants=paste("intfail at t=", t, sep=""), animals=NA, a_traitsM=NA, a_traitV=NA, p_traitsM=NA, p_traitV=NA)
      return(outputlist) 
      break()
    }
    
    comp <- (K-u %*% t(a_pops))/K
    ben <- (alpha %*% t(p_pops))/(alpha %*% t(p_pops)+r) #This is the problem param
    bento <- (alpha %*% t(p_pops)) #Total sum of all plant ineractions for each pollinator species. The problem is that we're not penalizing based on the environmental context. (We have set per capita interactions, but these are depending on the environmental context. Need to update it)
    attemptscaling <- alpha*t(t(alpha)/colSums(alpha)) #This is closer, but not quite right. It fits the alpha scaling based on the plant context, but doesn't actually take into account the abundance of competitors when looking at plant interactions
    scaled <- a_pops %*% attemptscaling
    rowSums(scaled)
    
    
    alpha[,1]*t(a_pops) #For plant species 1, how much benefit does it get from each pollinator species in each site? 
    
      #Find out dpop/dt
      a_change <- matrix(data=NA, nrow=nsites, ncol=n_animals)
      for(n in 1:n_animals){
        propOfplants <- 0
        nums <- a_pops[,n]*alpha[n,]*p_pops#/K#[,site] What if we standardize plant species to K? This seems like it might have worked to fix my problem
        for(i in 1:n_plants){
          #Numerator-plant reward given by each focal pollinator
          propOfplants <- propOfplants+nums[,i]/(r+rowSums(alpha[,i]*a_pops)) #Denominator-scaled by total available pollen
        }
        pollDeath <- mua
        #print(propOfplants)
        a_change[,n] <- lambda*((1-(o[,n]%*%t(a_pops))/K))*propOfplants# - pollDeath
        if(t %% 10==0){ #Trying to make the output object smaller, so only saving every 10th timestep
          propOutputtemp[[t %/% 10]][[n]] <- propOfplants #Save the benefit each pollinator receives from plants per capita as an entry in a list so we can compare expected vs realized specialism
        }
      } #Animal Growth Rate
      
      p_change <- matrix(data=NA, nrow=nsites, ncol=n_plants)
      for(n in 1:n_plants){
        plantCompetition <- (K-rowSums(t(u[n,]*t(p_pops))))/K#(K-sum(u[n,]*p_pops[n,]))/K
        #pollBenefits <-sum(alpha[,n]*a_pops*p_pops[n])/(r+sum(alpha[,n]*a_pops)) #NOTE: THIS IS MORE SIMILIAR TO BECKER. NOt sure why their p_pop is included here. 
        pollBenefits <-rowSums(t(alpha[,n]*t(a_pops)))/(r+rowSums(t(alpha[,n]*t(a_pops))))
        plantDeath <- mup
        p_change[,n] <- r*plantCompetition*pollBenefits #- plantDeath# + o*rnorm(n=1, mean=0, sd=0.1) 
        if(t %% 10==0){ #Trying to make the output object smaller, so only saving every 10th timestep
          pollBenefitsOutput[[t %/% 10]][[n]] <- pollBenefits #Save the benefit each pollinator receives from plants per capita as an entry in a list so we can compare expected vs realized specialism
        }
      }
      #Project to t+1
      for(i in 1:(n_plants*nsites)){
        #p_pops[i] <- sum(rpois(n=round(p_pops[i]),lambda = 1+p_change[i])) #Poisson birth process for plants
        p_pops[i] <- sum(rpois(n=round(p_pops[i]),lambda = exp(p_change[i])))
        #browser()
        #p_pops[i] <- temp
      }
      

      #p_pops <- p_pops*(1+p_change)
      #a_pops <- (a_change+1)*a_pops #Update Animal population
      for(i in 1:(n_animals*nsites)){
        #a_pops[i] <- sum(rpois(n=round(a_pops[i]),lambda = 1+a_change[i])) #Poisson birth process for animals
        new_a_pop <- try(sum(rpois(n=round(a_pops[i]), lambda = exp(a_change[i]))))
        if(class(new_a_pop)=='try-error'){
          lambda_fail <- TRUE
          next()
        }
        a_pops[i] <- new_a_pop
      }
      if(lambda_fail==TRUE){ #If we fail a lambda projection, we break our overall loop
        #outputlist <- list(warn=paste("lambdafail at t=", t, sep=""), plants=p_pops_output, animals=a_pops_output, a_traitsM=NA, a_traitV=NA, p_traitsM=NA, p_traitV=NA)
        outputlist <- list(warn=paste("lambdafail at t=", t, sep=""), plants=p_pops_output, animals=a_pops_output, a_traitsM=a_traitM, a_traitV=a_traitV, p_traitsM=p_traitM, p_traitV=p_traitV)
        return(outputlist) 
        
        break()
      }
      if(disptype=="positiveComp"){
        a_disprobs <- disprobmax*exp(a_change)/2
      }
      if(disptype=="negativeComp"){#Dispersal probability is inversely proportional to local growth rate
        a_disprobs <- 1-disprobmax*exp(a_change)/2
      }
      if(disptype=="neutralComp"){
        a_disprobs <- disprobmax/2 #Note: this creates a single scalar, rather than a matrix. Should work the same in the dispersal operation below, but just something to keep track of. 
      }
      
      a_disprobs[a_disprobs<0] <- 0 #Species that more than double have their dispersal prob set to 0
      a_disprobs[a_disprobs>disprobmax] <- disprobmax #Species that really do poorly can't have a dispersal probability above disprobmax
      #Dispersal (occurs after demographics for year)
      #Decide number of animal dispersers from each site of each species
      a_emms <- matrix(rbinom(n=length(a_pops), size=a_pops, prob=a_disprobs), nrow=nsites)
      a_pops <- a_pops-a_emms #Remove our immigrants from our pops for now so we don't count them redundantly
      
      #animals disperse across sites
      for(site in 1:nsites){
        for(spp in 1:n_animals){
          destin <- sample(factor_sites, size=a_emms[site, spp], replace=TRUE, prob=dexpsim[site,]/sum(dexpsim[site,])) #Sample our sites a number of times equal to the emmigrants from a particular site.
          a_pops[,spp] <- a_pops[,spp] + as.matrix(table(destin), ncol=1) #Table summarizes how many immigrants go to each site and then we add them to the population size. 
        }
      }
      
      p_disprobs <- disprobmax/2 #Plant dispersal probability is always a constant (can't realistically modify their behavior based on local conditions; changes in numerical response should already be reflected in growth rate)
      #Now our plants can disperse
      p_emms <- matrix(rbinom(n=length(p_pops), size=p_pops, prob=p_disprobs), nrow=nsites)
      p_pops <- p_pops-p_emms #Remove our immigrants from our pops for now so we don't count them redundantly
      
      #plants disperse across sites
      for(site in 1:nsites){
        for(spp in 1:n_plants){
          #browser()
          destin <- sample(factor_sites, size=p_emms[site, spp], replace=TRUE, prob=dexpsim[site,]/sum(dexpsim[site,])) #Sample our sites a number of times equal to the emmigrants from a particular site.
          p_pops[,spp] <- p_pops[,spp] + as.matrix(table(destin), ncol=1) #Table summarizes how many immigrants go to each site and then we add them to the population size. 
          #print(paste("site=", site, "; spp=", spp, sep=""))
        }
      }
      
      remove <- FALSE
      #Extinction?
      if(length(colSums(p_pops)[colSums(p_pops)<=e_thresh])>0){
        #print("global plant extinction")
        remove <- TRUE
      }
      if(remove==TRUE){
        #Plant Extinction
        #p_traitM <- p_traitM[which(colSums(p_pops)>e_thresh)] #Remove trait values from globally exinct plants
        #p_traitV <- p_traitV[which(colSums(p_pops)>e_thresh)]
        #p_pops <- p_pops[,which(colSums(p_pops)>e_thresh)] #Stop tracking populations of globally extinct plants
        #Animal Extinction
        #p_pops[,c(1:ncol(p_pops))[colSums(p_pops)<=e_thresh]] <- 0
        #remove <- FALSE
      }
      if(length(colSums(a_pops)[colSums(a_pops)<=e_thresh])>0){
        #print("global animal extinction")
        remove <- TRUE
        a_pops[,c(1:ncol(a_pops))[colSums(a_pops)<=e_thresh]] <- 0 #enforce extinction, rounding down any below 2.
      } 
      #Problem: keeping species in as 0's makes calculations run slower due to matrix issues. However, it makes it much easier to keep track of our total population sizes. Not sure which direction to solve this in. 
      #a_traitM <- a_traitM[which(colSums(a_pops)>e_thresh)] #Remove trait values from globally extinct plants
      #a_traitV <- a_traitV[which(colSums(a_pops)>e_thresh)]
      #a_pops <- a_pops[,which(colSums(a_pops)>e_thresh)] #Stop tracking populations of globally extinct plants
      
      #Do we invade? 
      if(rbinom(1, 1, invProb)==1){
        #print("invasion")
        invasion <- runif(1)
        invasion_vector <- c(rep(0,nsites))
        invasion_vector[sample(1:nsites, size=1)] <- invade_size #choose a random site to invade with the invade_size; this is set up to be the same invader population size for plants and polls, but this could easily be split up
        if(invasion >=0.5){
          #print("plant invades")
          p_pops <- cbind(p_pops, unlist(invasion_vector)) #new pop always invades a random site at invade_size; unlist is to get rid of that pesky name
          p_traitM <- c(p_traitM, runif(1)) #pick new trait mean from a uniform distribution [0,1]
          p_traitV <- c(p_traitV, runif(1, 0.05, 0.25)) #All set to be the same in Becker paper
        }else{
          #print("pollinator invades")
          a_pops <- cbind(a_pops, unlist(invasion_vector)) #new pop always invades a random site at invade_size
          a_traitM <- c(a_traitM, runif(1)) #pick new trait mean from a uniform distribution [0,1]
          a_traitV <- c(a_traitV, runif(1, 0.05, 0.25)) #Same for new trait variance
        }
      }
      n_animals <- ncol(a_pops)
      n_plants <- ncol(p_pops)
      
      if(n_animals < 1| n_plants < 1){
        #print("Network breakdown")
        break
      }
      #print(paste("polls=", n_animals, sep=""))
      #WmeanV <- rbind(WmeanV, sum((a_pops*a_traitV))/sum(a_pops))
      #meanV <- rbind(meanV, sum(a_traitV)/n_animals)
      #if(t %% 100==0)
      
      if((ncol(a_pops)+2)>ncol(a_pops_output)){ #change the output dimensions to account for a new species
        a_pops_output <- cbind(a_pops_output, c(rep(0,nrow(a_pops_output))))#add a 0 for every preceeding timestep
        colnames(a_pops_output)[ncol(a_pops_output)] <- colnames(a_pops)[ncol(a_pops)] # rename based on the newest animal species
      }
      a_pops_temp <- cbind(c(1:nsites), rep(t,nsites),a_pops)
      a_pops_output <- rbind(a_pops_output, a_pops_temp)
      
      
      if((ncol(p_pops)+2)>ncol(p_pops_output)){ #change the output dimensions to account for a new species
        p_pops_output <- cbind(p_pops_output, c(rep(0,nrow(p_pops_output))))
        colnames(p_pops_output)[ncol(p_pops_output)] <- colnames(p_pops)[ncol(p_pops)] #  rename based on the newest plant species
      }
      p_pops_temp <- cbind(c(1:nsites), rep(t,nsites),p_pops)
      p_pops_output <- rbind(p_pops_output, p_pops_temp)
      
      #pops <- rbind(pops, pops_temp)
      if(t %% 500==0)
        print(t)
      
      P_WTraitMean_site <- data.frame(site=1:nsites, t=t, WVMean=NA)
      A_WTraitMean_site <- data.frame(site=1:nsites, t=t, WVMean=NA)
      for(i in 1:nsites){
        P_WTraitMean_site$WVMean[i] <- sum(p_pops[i,]*p_traitV)/sum(p_pops[i,])
        A_WTraitMean_site[i] <- sum(a_pops[i,]*a_traitV)/sum(a_pops[i,])
      }
      
      P_WTraitMean <- sum(colSums(p_pops)*p_traitV/sum(p_pops))
      A_WTraitMean <- sum(colSums(a_pops)*a_traitV/sum(a_pops))
      
      P_WTraitMean_output <- rbind(P_WTraitMean_output, P_WTraitMean_site)
      A_WTraitMean_output <- rbind(A_WTraitMean_output, P_WTraitMean_site)
      
      P_TraitVTot <- rbind(P_TraitVTot, P_WTraitMean)
      A_TraitVTot <- rbind(A_TraitVTot, A_WTraitMean)
  }
  outputlist <- list(plants=p_pops_output, animals=a_pops_output, a_traitsM=a_traitM, a_traitV=a_traitV, p_traitsM=p_traitM, p_traitV=p_traitV, RSpecPoll=propOutputtemp, RSpecPlant=pollBenefitsOutput)
  return(outputlist) #Specify return object
}

