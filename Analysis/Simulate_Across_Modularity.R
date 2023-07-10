
library(scales)
library(parallel)
library(igraph)

library(igraph)
nsites <- 100


A <- matrix(data=runif(nsites^2, 0,1),ncol=nsites)
Asym <- A %*% t(A)
class(Asym/max(Asym))

library(igraph)
nsites <- 100



wi <-0.2 #Bernoulli probability of link within module
ac <-0.05 #Bernoulli probability of link between module
pref <- matrix(data=c(wi, ac, ac, ac, ac, #Create a preference matrix from the above vals
             ac, wi, ac, ac, ac,
             ac, ac, wi, ac, ac,
             ac, ac, ac, wi, ac,
             ac, ac, ac, ac, wi), ncol=5, byrow = TRUE)
a <- igraph::sample_sbm(nsites, pref.matrix=pref, block.sizes = rep(20, nsites/20)) #Use sample_sbm to make a random graph

com <- igraph::fastgreedy.community(a) #Find community identity

igraph::modularity(a, membership=com$membership, weights=E(a)$weights) #Get resulting modularity
hist(igraph::degree(a))
igraph::vertex.connectivity(a)

coords <- layout_nicely(a)
eucdist <- dist(coords, diag=T, upper=T)
dexpdist <- dexp(eucdist, rate=1)
dexpsim <- 1/(1+as.matrix(dexpdist)) #This has a full distance matrix. Do I only assign the values for the sites connected in my original graph, or do I let them all have pairwise connections?

plot(a, layout=coords)
E(a)$weights[1]

val <- E(a)[1]

dists <- distance_table(a, directed = FALSE)
dists <- igraph::distances(a)
hist(dists)



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

save(list_results, file="simulationNegative_30sites.Rda")
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

save(list_results, file="simulationPositive_30sites.Rda")
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

save(list_results, file="simulationNeutral_30sites.Rda")
print("Done with neutral")

save(time1, time2, time3, file="timing.Rda") #Save out output time objects
