
library(scales)
library(parallel)
library(igraph)

library(igraph)
library(Matrix)
library(simplifyNet)
nsites <- 100


A <- matrix(data=runif(nsites^2, 0,1),ncol=nsites)
Asym <- A %*% t(A)
class(Asym/max(Asym))

nsites <- 40



wi <-0.4 #Bernoulli probability of link within module
ac <-0.01 #Bernoulli probability of link between module
pref <- matrix(data=rep(ac, (nsites/20)^2), ncol=nsites/20, byrow = TRUE)
diag(pref) <- wi
a <- igraph::sample_sbm(nsites, pref.matrix=pref, block.sizes = rep(20, nsites/20)) #Use sample_sbm to make a random graph
com <- igraph::fastgreedy.community(a) #Find community identity
V(a)$color <- com$membership

coords <-layout_nicely(a)
plot(a, layout=coords)


empty <- matrix(data=1, nrow=nsites, ncol=nsites)
diag(empty) <- 0
graph <- graph_from_adjacency_matrix(adjmatrix = empty, mode=c("undirected"))
E(graph)$weights <- dexpsim[t(combn(nrow(dexpsim), 2))]
V(graph)$color <- com$membership
igraph::fastgreedy.community(graph, weights=E(graph)$weights)
igraph::fastgreedy.community(g2, weights=E(g2)$weights)

g2 <- delete.edges(graph, which(E(graph)$weights > quantile(E(graph)$weights,0.2))) #This should be trimming down to only the 20% of closest links, but seems random. Looks like link weight is not being assigned correctly indexed.


plot(g2, layout=coords)
plot(a, layout=coords)
for(i in 1:nsites){
  if(all(V(a)$color[neighbors(a, v=i)]==V(a)$color[i])==FALSE){ #Check if this node has connections outside its module

    E(a)[V(a)[color==1] %--% V(a)[color==2]]
    
    E(graph)$weights[E(graph)==E(a)[V(a)[color==1] %--% V(a)[color==2]]] #extract the weights of the intermodules
  }
}
as_adjacency_matrix(graph)
E(graph)$weight <- E(graph)$weights

S <- simplifyNet::bestpath(graph)

sg = simplifyNet::net.as(S, net.to="igraph", directed=FALSE)
plot(sg)
plot(graph)
igraph::ecount(sg)/igraph::ecount(graph)#fraction of edges in the sparsifier


#Generate random ER graph with uniformly random edge weights
g = igraph::erdos.renyi.game(50, 0.1)
igraph::E(g)$weight <- runif(length(igraph::E(g)))
#Sparsify g via bestpath
S = simplifyNet::bestpath(g, directed = FALSE, associative = TRUE) #Show edge list conversion
sg = simplifyNet::net.as(S, net.to="igraph", directed=FALSE)
igraph::ecount(sg)/igraph::ecount(graph)#fraction of edges in the sparsifier using bestpath


effR <- simplifyNet::EffR(graph)
Eff <- EffRSparse(graph, 100, effR, 24601)
Effg <- simplifyNet::net.as(Eff, net.to="igraph", directed=F)
igraph::ecount(Effg)/igraph::ecount(graph)
modules <- igraph::fastgreedy.community(Effg, weights=E(Effg)$weights) #Similiar to my problem before, the new sparsified graph shows different modules than the sbm graph that it was buil upon :/

V(Effg)$color <- modules$membership

modularity(Effg, weights=E(Effg)$weights, membership=modules$membership) 

plot(Effg)

EffRSparse(network, q, effR, 24601, n)
get.adjacency.sparse(graph)
E(graph)$weights

igraph::modularity(g2, membership=com$membership, weights=E(g2)$weights) #Get resulting modularity
hist(igraph::degree(a))
igraph::vertex.connectivity(a)
plot(a)

#####################################################################################



coords <- layout_nicely(a)
#coords <- igraph::layout_with_kk(a)
coords
coords[,1] <- coords[,1]/max(abs(coords[,1])) #rescale to between -1 to 1 so dexp is comparable to my random graphs
coords[,2] <- coords[,2]/max(abs(coords[,2]))#rescale to between -1 to 1 so dexp is comparable to my random graphs
eucdist <- dist(coords, diag=T, upper=T)#r
dexpdist <- dexp(eucdist, rate=2.5) #This rate parameter is somewhat arbitrarily chosen. Used 10 in other simulations, but relized this gives a bunch of nearly 1 connectivities.
dexpsim<- 1/(1+as.matrix(dexpdist)) #This has a full distance matrix. Do I only assign the values for the sites connected in my original graph, or do I let them all have pairwise connections?
hist(x)
hist(x[x>quantile(0.8, x)])



empty <- matrix(data=1, nrow=nsites, ncol=nsites)
diag(empty) <- 0
dexpsim[upper.tri(dexpsim)==TRUE]

graph <- graph_from_adjacency_matrix(adjmatrix = empty, mode=c("undirected"))




triu(dexpsim)[triu(dexpsim)!=1][1:3]
dexpsim[upper.tri(dexpsim)==TRUE][1:3]
dexpsim[1:4,1:4]
?triu
result <- dexpsim[t(combn(nrow(dexpsim), 2))] #This code extracts the non-diagonal upper triangle values of the distance matrix rowwise (so they're in an order that can be assigned to the rowweights)

tail(a)
E(graph)$weights <- dexpsim[t(combn(nrow(dexpsim), 2))] #This code extracts the non-diagonal upper triangle values of the distance matrix rowwise (so they're in an order that can be assigned to the rowweights)

fastgreedy.community(graph, weights=E(graph)$weights)
cluster_edge_betweenness(graph, weights=E(graph)$weights)

#Problem: We create a modular graph that contains 5 submodules using sample_sbm(). However, this only creates 0-1 edges (no weight)
#To make weights approximating a real-world landscape, we layout_nicely this graph, turn the coordinates into a euclidean distance matrix, run that through a exponential decay
#And now we have a new graph, fully connected, where edgeweights are equivalent to the connectivity between sites.
#The only downside to this is that due to this process, we no longer really have the same number of modules, and I'm not sure how much variability in modularity we have to control
#When we degrade networks from a high modularity, the step probably needs done at the last step. I guess I'll do this by randomly swapping network link weights?
quantile(E(graph)$weights,0.8)
g2 <- delete.edges(graph, which(E(graph)$weights < quantile(E(graph)$weights,0.95))) #This should be trimming down to only the 20% of closest links, but seems random. Looks like link weight is not being assigned correctly indexed.
hist(E(graph)$weights)
plot(a, layout=coords)
plot(g2, layout=coords)
igraph::modularity(graph, membership=fastgreedy.community(graph, weights=E(graph)$weights)$membership, weights=E(graph)$weights) #Get resulting modularity. 

plot(graph, layout=coords)
E(a)$weights[1]

val <- E(a)[1]

dists <- distance_table(a, directed = FALSE)
dists <- igraph::distances(a)

range(dexpsim)

factor_sites <- as.factor(1:nsites) #This is important for the dispersal function down the road
coords <- data.frame(x=runif(nsites), y=runif(nsites))
eucdist <- dist(coords, diag=T, upper=T)
dexpdist <- dexp(eucdist, rate=20)
dexpsim <- 1/(1+as.matrix(dexpdist)) #Converting our distance matrix to a similarity score (easier for me to then turn into a pmf).
range(dexpsim)
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
