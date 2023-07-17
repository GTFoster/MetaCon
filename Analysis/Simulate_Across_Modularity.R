
library(scales)
library(parallel)
library(igraph)

library(igraph)
library(Matrix)
library(simplifyNet)
nsites <- 100


wi <-0.4 #Bernoulli probability of link within module
ac <-0.01 #Bernoulli probability of link between module
pref <- matrix(data=rep(ac, (nsites/25)^2), ncol=nsites/25, byrow = TRUE) #Create sbm preference matrix according to above probabilities
diag(pref) <- wi 
a <- igraph::sample_sbm(nsites, pref.matrix=pref, block.sizes = rep(25, nsites/25)) #Use sample_sbm to make a random graph
com <- igraph::fastgreedy.community(a) #Find community identity
V(a)$color <- com$membership #assign it to color

coords <- layout_nicely(a) #use layout-nicely to enforce 2D spatial orientation
coords[,1] <- coords[,1]/max(abs(coords[,1])) #rescale to between -1 to 1 so dexp is comparable to my random graphs
coords[,2] <- coords[,2]/max(abs(coords[,2]))#rescale to between -1 to 1 so dexp is comparable to my random graphs
eucdist <- dist(coords, diag=T, upper=T) #Create distance matrix from coordinates. 
dexpdist <- dexp(eucdist, rate=2.5) #This rate parameter is somewhat arbitrarily chosen. Used 10 in other simulations, but relized this gives a bunch of nearly 1 connectivities.
dexpsim <- as.matrix(dexpdist/max(dexpdist)) #NOTE: Changed from previous formulation. Realized that a distance score wasn't great to combine with network weighting below
#plot(a, layout=coords)

empty <- matrix(data=1, nrow=nsites, ncol=nsites) #Create fully-connected adjacency matrix
diag(empty) <- 0 #remove self-loops
graph <- graph_from_adjacency_matrix(adjmatrix = empty, mode=c("undirected")) #turn it into a graph
E(graph)$weight <- dexpsim[t(combn(nrow(dexpsim), 2))] #assign weights (connectivitied) based on the exponential distance kernel

V(graph)$color <- com$membership #Assign membership using the ORIGINAL fastgreedy optimization for visualization

#plot(graph, layout=coords) #This graph is super dense; going to reduce using bestpath sparsification (package by Kramer, 2022, method by Spielman et al., 2013)
#g2 <- delete.edges(graph, which(E(graph)$weight < quantile(E(graph)$weight, .8))) #Trim down to only the closest 20% links-old, really coarse sparsification method

#plot(a, layout=coords)

s <- simplifyNet::bestpath(graph, associative = TRUE) #Conduct sparsification
sg <- simplifyNet::net.as(s, net.to="igraph", directed=FALSE) #return sparsified graph as igraph object
#plot(sg, layout=coords)
igraph::ecount(sg)/igraph::ecount(graph)#fraction of edges in the sparsified network

classifysg <- fastgreedy.community(sg, weights=E(sg)$weight) #Recalculate community membership. Since we've changed our structure by fully connecting and then sparsifying, we need to reassign communities to avoid artifically reducing modularity score
modularity(sg, membership = classifysg$membership) #Quite modular network. Next step is to degrade its modularity by randomly reassigning links


sgdeg <- sg 
modularity(sg, membership = classifysg$membership) #What's our starting modularity?
output <- NULL
ndegrades <- 1000 #How many links should we redistribute?
failcount <- 0
i <- 1
brake <- 0

while(i<500){
  brake <- brake +1 #Make sure we don't get stuck here
  if(brake > ndegrades*1000){
    print("too many failed reassignments-stuck in while loop")
    break()
  }
  new_edge <- sample(1:length(V(sgdeg)), size=2, replace=FALSE) #Randomly choose two new vertices to connect
  if(are.connected(sgdeg, v1=new_edge[1], v2=new_edge[2])==TRUE){ #Check if that new edge already exists in the network
    next() #If it does already exist, move to the next iteration. NEED TO MAKE THIS TRY AGAIN RATHER THAN JUST GIVE UP
    }
  
  old_edge <- sample(1:length(E(sgdeg)), size=1) #Choose an old edge to delete
  old_weight <- E(sgdeg)$weight[old_edge] #extract its weight
  sgdeg <- delete_edges(sgdeg, edges=old_edge) #Delete it from the graph
  
  sgdeg <- add_edges(sgdeg, edges=new_edge, weight=old_weight) #Add in the new edge, with its new weight equal to the deleted (reassigning existing edge set)
  temp <- data.frame(t=i,mode=modularity(sgdeg, membership=classifysg$membership)) #return the number of edges reassigned, plus the new modularity score
  output <- rbind(output, temp) 
  i <- i+1
}

plot(x=output$t, y=output$mode)
#####################################################################################



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
