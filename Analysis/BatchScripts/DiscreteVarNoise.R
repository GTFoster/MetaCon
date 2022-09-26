output <- NULL
for(q in 1:200){
  n_plants <- 3
  n_animals <- 3
  
  p_pops <- c(.10, .05, 0.25)
  p_traitM <- c(0.25, 0.99, 0.5)
  p_traitV <- c(0.1, 0.1, 0.1) #All set to be the same in Becker paper
  
  #a_pops[1:n_animals] <- c(5, 10, 20)
  a_pops <- c(.5, .1, .2)
  a_traitM <- c(0.1, 0.66, 0.8)
  a_traitV <- c(0.25, 0.5, 0.1)
  
  r <- 1.0
  mup <- .2
  mua <- 0.1
  if(q<101){ #100 runs of each 0 value
    o <- 0.1
  }
  if(q>100){
    o <- 1
  }

  num_timeSteps <- 400
  p_rich <- NULL
  a_rich <- NULL
  WmeanV <- NULL
  meanV <- NULL
  
  for(t in 1:num_timeSteps){
    p_rich <- rbind(p_rich, n_plants)
    a_rich <- rbind(a_rich, n_animals)
    #Check if community is same as in t-1
    #if so; recalculate A and u
    alpha <- matrix(nrow=n_animals, ncol=n_plants) #Create an empty matrix for plant-pollinator partnerships
    for(ia in 1:n_animals){
      for(ip in 1:n_plants){
        alpha[ia, ip] <- integrate(min.f1f2, -Inf, Inf, mu1=p_traitM[ip], mu2=a_traitM[ia], sd1=p_traitV[ip], sd2=a_traitV[ia])$value #Integrate our minimum equation over all numbers; output is the total area overlapping both curves, with max value of perfectly overlapping curves as 1
      }
    }
    alpha[alpha<0.05] <- 0
    
    u <- matrix(data=NA, nrow=n_plants, ncol=n_plants)
    
    for(p1 in 1:n_plants){
      for(p2 in 1:n_plants){
        u[p1, p2] <- integrate(min.f1f2, -Inf, Inf, mu1=p_traitM[p1], mu2=p_traitM[p2], sd1=p_traitV[p1], sd2=p_traitV[p2])$value #Integrate our minimum equation over all numbers; output is the total area overlapping both curves, with max value of perfectly overlapping curves as 1
      }
    }
    #Find out dpop/dt
    a_change <- vector(mode="numeric", length=n_animals) #Plant growth rate
    for(n in 1:n_animals){
      propOfplants <- 0
      nums <- alpha[n,]*p_pops
      for(i in 1:n_plants){
        #Numerator-plant reward given by each focal pollinator
        propOfplants <- propOfplants+nums[i]/(r+sum(alpha[,i]*a_pops)) #Denominator-scaled by total available pollen
      }
      pollDeath <- mua*a_pops[n]
      a_change[n] <- lambda*propOfplants - pollDeath + o*rnorm(n=1, mean=0, sd=0.1)+o*rnorm(n=1, mean=0, sd=0.1)
    } #Animal Growth Rate
    
    p_change <- vector(mode="numeric", length=n_plants) #Plant growth rate
    for(n in 1:n_plants){
      plantCompetition <- 1-sum(u[n,]*p_pops[n])
      pollBenefits <-sum(alpha[,n]*a_pops*p_pops[n])/(r+sum(alpha[,n]*a_pops))
      plantDeath <- mup*p_pops[n]
      
      p_change[n] <- r*plantCompetition*pollBenefits - plantDeath + o*rnorm(n=1, mean=0, sd=0.1) 
    }
    #Project to t+1
    p_pops <- p_pops*(1+p_change)
    a_pops <- (a_change+1)*a_pops #Update Animal population
    #Extinction?
    if(length(p_pops[p_pops<0.05])>0){
      print("plant extinction")
    }
    #Plant Extinction
    p_pops <- p_pops[which(p_pops>0.05)]
    p_traitM <- p_traitM[which(p_pops>0.05)]
    p_traitV <- p_traitV[which(p_pops>0.05)]
    #Animal Extinction
    if(length(a_pops[a_pops<0.05])>0){
      print("animal extinction")
    }
    a_pops <- a_pops[which(a_pops>0.05)]
    a_traitM <- a_traitM[which(a_pops>0.05)]
    a_traitV <- a_traitV[which(a_pops>0.05)]  
    #Do we invade? 
    if(rbinom(1, 1, 0.05)==1){
      #print("invasion")
      invasion <- runif(1)
      if(invasion >=0.5){
        print("plant invades")
        p_pops <- c(p_pops, 0.1) #new pop always invades at 0.1
        p_traitM <- c(p_traitM, runif(1)) #pick new trait mean from a uniform distribution [0,1]
        p_traitV <- c(p_traitV, 0.1) #All set to be the same in Becker paper
      }else{
        print("pollinator invades")
        a_pops <- c(a_pops, 0.1) #new pop always invades at 0.1
        a_traitM <- c(a_traitM, runif(1)) #pick new trait mean from a uniform distribution [0,1]
        a_traitV <- c(a_traitV, runif(1)) #Same for new trait variance
      }
    }
    n_animals <- length(a_pops)
    n_plants <- length(p_pops)
    if(n_animals < 1| n_plants < 1){
      break
    }
    
    WmeanV <- rbind(WmeanV, sum((a_pops*a_traitV))/n_animals)
    meanV <- rbind(meanV, sum(a_traitV)/n_animals)
  }
  
  temp <- data.frame(Wv=WmeanV, mV=meanV, run=1, Arich=a_rich, Prich=p_rich, time=1:nrow(WmeanV), O=o)
  output <- rbind(output, temp)
}
save(output, file="DiscreteVarNoise.Rda")
