
library(scales)
library(parallel)




nsites <- 5
factor_sites <- as.factor(1:nsites) #This is important for the dispersal function down the road
coords <- data.frame(x=runif(nsites), y=runif(nsites))
eucdist <- dist(coords, diag=T, upper=T)
dexpdist <- dexp(eucdist, rate=10)
dexpsim <- 1/(1+as.matrix(dexpdist)) #Converting our distance matrix to a similarity score (easier for me to then turn into a pmf).
#Found this particular metric from here: https://stats.stackexchange.com/questions/158279/how-i-can-convert-distance-euclidean-to-similarity-score





source(file="./DispersalSimulationAnComp.R") #load in our dispersal simulation function. rOxygen-style description available in file.




#list_results <- list()
#for(i in 1:5){
#list_results[[i]] <- try(runDispersalSim(nsites=nsites, disptype="negativeComp",n_plants=5, n_animals=5, dexpsim=dexpsim, r=0.5, mup=0.1, mua=0.1, o=0.1, lambda=0.9, K=100, e_thresh = 2, invade_size = 5, disprob = 0.2, num_timeSteps = 3000, invProb=0.5))
#}





tictoc::tic()
num_iterations <- 100
#num_iterations <- 4

# Set up the cluster for parallelization
cl <- makeCluster(detectCores()-2, outfile="")
              
# Define a function to run the simulation
runSimulation <- function() {
  x <- try(runDispersalSim(nsites = nsites, disptype = "negativeComp", n_plants = 5, n_animals = 5, dexpsim = dexpsim, r = 0.5, mup = 0.1, mua = 0.1, o = 0.1, lambda = 0.9, K = 500, e_thresh = 2, invade_size = 5, disprobmax = 0.2, num_timeSteps = 3000, invProb = 0.05))
  return(x)
  }

clusterExport(cl, c("nsites", "runSimulation", "runDispersalSim", "dexpsim"))

# Run the simulations in parallel
list_results <- clusterApplyLB(cl, 1:num_iterations, function(i) {
  runSimulation()
})

# Stop the cluster
stopCluster(cl)
time1 <- tictoc::toc()

save(list_results, file="simulationNegative.Rda")
print("Done with negative")



tictoc::tic()
cl <- makeCluster(detectCores()-2, outfile="")
runSimulation <- function() {
  try(runDispersalSim(nsites = nsites, disptype = "positiveComp", n_plants = 5, n_animals = 5, dexpsim = dexpsim, r = 0.5, mup = 0.1, mua = 0.1, o = 0.1, lambda = 0.9, K = 500, e_thresh = 2, invade_size = 5, disprob = 0.2, num_timeSteps = 3000, invProb = 0.05))
}
clusterExport(cl, c("nsites", "runSimulation", "runDispersalSim", "dexpsim"))

# Run the simulations in parallel
list_results <- clusterApplyLB(cl, 1:num_iterations, function(i) {
  try(runSimulation())
})

stopCluster(cl)
time2 <- tictoc::toc()

save(list_results, file="simulationPositive.Rda")
print("Done with positive")



tictoc::tic()
cl <- makeCluster(detectCores()-2, outfile="")
runSimulation <- function() {
  try(runDispersalSim(nsites = nsites, disptype = "neutralComp", n_plants = 5, n_animals = 5, dexpsim = dexpsim, r = 0.5, mup = 0.1, mua = 0.1, o = 0.1, lambda = 0.9, K = 500, e_thresh = 2, invade_size = 5, disprob = 0.2, num_timeSteps = 3000, invProb = 0.05))
}
clusterExport(cl, c("nsites", "runSimulation", "runDispersalSim", "dexpsim"))

# Run the simulations in parallel
list_results <- clusterApplyLB(cl, 1:num_iterations, function(i) {
  runSimulation()
})

runSimulation()

stopCluster(cl)
time3 <- tictoc::toc()

save(list_results, file="simulationNeutral.Rda")
print("Done with neutral")


