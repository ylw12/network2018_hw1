#################################################################
##################### Network2018 Homework1 #####################
#####################       Luwei Ying      #####################
#################################################################

rm(list = ls())
library('igraph')
library('sna')
library('statnet')

setwd('~/Dropbox/Summer_Courses/workshopWashU/network2018_hw1')

# -------- 1. Nigeria Data Processing
load("nigeria.rda")
head(nigeria)
range(nigeria$year)

# Turn the data into a matrix
years <- sort(unique(nigeria$year))
countries <- with(nigeria, intersect(sender, receiver))
n <- length(countries)
adjMat <- matrix(0, nrow = n, ncol = n, dimnames = list(countries, countries))

conflictlist  <- lapply(years, function(t){
  slice <- nigeria[nigeria$year == t,]
  CflCases <- slice[slice$conflict == 1,]
  for(i in 1:nrow(CflCases)){
    adjMat[as.character(CflCases[i,]$sender), as.character(CflCases[i,]$receiver)] <- 1
  }
  return(adjMat)
})
names(conflictlist) <- years

# -------- 2. Measurements & Community detection

# Create a graph object from the adjacency matrix
graph <- graph_from_adjacency_matrix(conflictlist[[1]], 
                                     mode='directed', 
                                     weighted=TRUE,
                                     diag=FALSE
)

# Influential actor (2000)
which.max(degree(graph))
# Using degree, the police of Nigeria is the most influencial actor. 
# There are 21 edges connected to the police of Nigeria. 
# That means it is the actor being involved in the most conflicts.

which.max(eigen_centrality(graph, directed = TRUE)$vector)
# Using eigenvector centrality, Fulani Militia is the most influencial actor.
# The high eigenvector score of Fulani Militia means it is connected to many 
# other military organizations who themselves have high scores.

# Create clusters based on structural equivalence
eq <- equiv.clust(conflictlist)
plot(eq,hang=-1)
dev.off()

# Try k = 2
Conflictblocks2 <- blockmodel(conflictlist, eq, k = 2)
NodeCls2 <- Conflictblocks2$block.membership
countries[c(17, 21)] # Military (Nigeria), Police (Nigeria)

# Try k = 3
Conflictblocks3 <- blockmodel(conflictlist, eq, k = 3)
NodeCls3 <- Conflictblocks3$block.membership
countries[17] # Military (Nigeria)
countries[21] # Police (Nigeria)

# Try k = 4
Conflictblocks4 <- blockmodel(conflictlist, eq, k = 4)
NodeCls4 <- Conflictblocks4$block.membership
countries[6] # Fulani Militia
countries[17] # Military (Nigeria)
countries[21] # Police (Nigeria)

# Try k = 5
Conflictblocks5 <- blockmodel(conflictlist, eq, k = 5)
NodeCls5 <- Conflictblocks5$block.membership
countries[6] # Fulani Militia
countries[17] # Military (Nigeria)
countries[21] # Police (Nigeria)
countries[37] # Boko Haram

# Based on the cluster dendrogram, there seems no need to try a larger k
# once the two government institutions are detected.
# We'll cross validate the above 4 models.




# -------- 3. ERGMs

# To do the cross-sectional analysis, look at only 2000
Cfl2000 <- conflictlist[[1]]
diag(Cfl2000) <- NA
CflNet <- as.network.matrix(Cfl2000)

# Hypothesis: Reciprocity affects the network.
Mod <- ergm(CflNet ~ edges + mutual)
summary(Mod)

# Discuss: "mutual" is positive and significant in the model.
# That means when one military orgnization in nigeria attacks the 
# other, the other one is likely to fight back.

# check convergence
mcmc.diagnostics(Mod)

