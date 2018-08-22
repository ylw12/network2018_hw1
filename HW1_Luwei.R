#################################################################
##################### Network2018 Homework1 #####################
#####################       Luwei Ying      #####################
#################################################################

rm(list = ls())
library('cvTools')
library('igraph')
library('network')
library('PRROC')
library('sna')
library('statnet')
library('viridis')

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

g1 <- network(conflictlist[[1]], directed = T)

degree <- sna::degree(g1, cmode = 'outdegree')
degree <- ifelse(degree != 0, 1, 0)

ClustNum <- as.data.frame(cbind(NodeCls2, NodeCls3, NodeCls4, NodeCls5, degree))

# Set-up cross validation
folds <- cvFolds(nrow(ClustNum), K = 10)
ClustNum$Pred2 <- rep(0, nrow(ClustNum))
ClustNum$Pred3 <- rep(0, nrow(ClustNum))
ClustNum$Pred4 <- rep(0, nrow(ClustNum))
ClustNum$Pred5 <- rep(0, nrow(ClustNum))

set.seed(1)
for (i in 1:10){
  train <- ClustNum[folds$subsets[folds$which != i], ]
  validation <- ClustNum[folds$subsets[folds$which == i], ]
  
  glm2 <- glm(degree ~ NodeCls2, data = train, family = binomial(link = 'logit'))
  pred2 <- predict(glm2, newdata = validation, type = 'response')
  ClustNum[folds$subsets[folds$which == i], ]$Pred2 <- pred2
  
  glm3 <- glm(degree ~ NodeCls3, data = train, family = binomial(link = 'logit'))
  pred3 <- predict(glm3, newdata = validation, type = 'response')
  ClustNum[folds$subsets[folds$which == i], ]$Pred3 <- pred3
  
  glm4 <- glm(degree ~ NodeCls4, data = train, family = binomial(link = 'logit'))
  pred4 <- predict(glm4, newdata = validation, type = 'response')
  ClustNum[folds$subsets[folds$which == i], ]$Pred4 <- pred4
  
  glm5 <- glm(degree ~ NodeCls5, data = train, family = binomial(link = 'logit'))
  pred5 <- predict(glm5, newdata = validation, type = 'response')
  ClustNum[folds$subsets[folds$which == i], ]$Pred5 <- pred5
}

for (i in 6:9){
  FG <- ClustNum[ClustNum$degree == 1, i]
  BG <- ClustNum[ClustNum$degree == 0, i]
  ROC[[i-5]] <- roc.curve(FG, BG, curve = T)
  PR[[i-5]] <- pr.curve(FG, BG, curve = T)
}
c(ROC[[1]]$auc, ROC[[2]]$auc, ROC[[3]]$auc, ROC[[4]]$auc)
c(PR[[1]]$auc.integral, PR[[2]]$auc.integral, PR[[3]]$auc.integral, PR[[4]]$auc.integral)


# Based on the results, choose k = 2 to visualize.
tiesSum = apply(graph[], 1, sum)
# condition size based on # of ties
V(graph)$size <- sqrt(tiesSum + 1)

# color
colnames(graph[])

bcols <- c(rep("blue", 16), "yellow", rep("blue", 3), "yellow", rep("blue", 16))

# plot!
par(mar=c(0,0,0,0))
plot(graph,
     vertex.label=V(graph)$name, 
     vertex.size=5*V(graph)$size,
     vertex.color = bcols, # change color of nodes
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .75, # change size of labels to 75% of original size
     edge.arrow.size = .25,
     edge.curved=.25, # add a 25% curve to the edges
     edge.color="grey60", # change edge color to grey
     layout=layout_on_grid
)


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

