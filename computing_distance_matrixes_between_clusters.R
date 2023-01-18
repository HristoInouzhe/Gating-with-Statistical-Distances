require(cubature)
require(transport)
require(kernlab)
require(flowMap)

source("auxiliary_functions.R")

# X and Y are defined in plots_figure_barycenter.R

data.X <- data.frame(X, id = X.clustered$cluster)
data.Y <- data.frame(Y, id = Y.clustered$cluster)

X.clust.labels <- names(table(data.X$id))
Y.clust.labels <- names(table(data.Y$id))

XYClust.Dist <- array(0, dim = c(length(X.clust.labels), length(Y.clust.labels)))
for(i in X.clust.labels){
  i.indexes <- which(data.X$id == i)
  i.n <- length(i.indexes)
  i.prob <- rep(1 / i.n, i.n)
  for(j in Y.clust.labels){
    j.indexes <- which(data.Y$id == j)
    j.n <- length(j.indexes)
    j.prob <- rep(1 / j.n, j.n)
    D2ij <- (as.matrix(dist(rbind(data.X[i.indexes, 1:2], data.Y[j.indexes, 1:2])))^2)[1:i.n, (i.n + 1):(i.n + j.n)]
    ij.ot.solution <- transport::transport(i.prob, j.prob, D2ij)  
    amount.mass <- array(0, dim = c(i.n, j.n))
    amount.mass[as.matrix(ij.ot.solution[, 1:2])] <- ij.ot.solution$mass
    w2dist <- 0
    for(k in 1:i.n){
      w2dist <- w2dist + D2ij[k, ] %*% amount.mass[k, ] 
    }
    XYClust.Dist[as.numeric(i), as.numeric(j)] <- sqrt(w2dist)
  }
}
XYClust.Dist.Norm <- XYClust.Dist / max(XYClust.Dist)

XYClust.Dist.MMD.vanilla <- array(0, dim = c(length(X.clust.labels), length(Y.clust.labels)))
k.kernel <- "vanilladot"
k.kpar <- list()
for(i in X.clust.labels){
  i.indexes <- which(data.X$id == i)
  for(j in Y.clust.labels){
    j.indexes <- which(data.Y$id == j)
    XYClust.Dist.MMD.vanilla[as.numeric(i), as.numeric(j)] <- kernlab::kmmd(x = as.matrix(data.X[i.indexes, 1:2]),
                                                                    y = as.matrix(data.Y[j.indexes, 1:2]),
                                                                    kernel = k.kernel,
                                                                    kpar = k.kpar, ntimes = 1)@mmdstats[2]
  }
}
XYClust.Dist.MMD.vanilla.Norm <- XYClust.Dist.MMD.vanilla / max(XYClust.Dist.MMD.vanilla)

XYClust.Dist.MMD.gauss <- array(0, dim = c(length(X.clust.labels), length(Y.clust.labels)))
k.kernel <- "rbfdot"
k.kpar <- list(0.001) 
for(i in X.clust.labels){
  i.indexes <- which(data.X$id == i)
  for(j in Y.clust.labels){
    j.indexes <- which(data.Y$id == j)
    XYClust.Dist.MMD.gauss[as.numeric(i), as.numeric(j)] <- kernlab::kmmd(x = as.matrix(data.X[i.indexes, 1:2]),
                                                                    y = as.matrix(data.Y[j.indexes, 1:2]),
                                                                    kernel = k.kernel,
                                                                    kpar = k.kpar, ntimes = 1)@mmdstats[2]
  }
}
XYClust.Dist.MMD.gauss.Norm <- abs(XYClust.Dist.MMD.gauss) / max(abs(XYClust.Dist.MMD.gauss))

XYClust.Dist.FR <- array(0, dim = c(length(X.clust.labels), length(Y.clust.labels)))
for(i in X.clust.labels){
  i.indexes <- which(data.X$id == i)
  for(j in Y.clust.labels){
    j.indexes <- which(data.Y$id == j)
    XYClust.Dist.FR[as.numeric(i), as.numeric(j)] <- flowMap::getFR(data.X[i.indexes, 1:2],
                                                                    data.Y[j.indexes, 1:2])$ww
  }
}
XYClust.Dist.FR <- abs(XYClust.Dist.FR)
XYClust.Dist.FR.Norm <- XYClust.Dist.FR / max(XYClust.Dist.FR)

XYClust.Dist.Hel <- array(0, dim = c(length(X.clust.labels), length(Y.clust.labels)))
for(i in X.clust.labels){
  i.indexes <- which(data.X$id == i)
  density.i <- ks::kde(data.X[i.indexes, 1:2], density = TRUE) 
  for(j in Y.clust.labels){
    j.indexes <- which(data.Y$id == j)
    density.j <- ks::kde(data.Y[j.indexes, 1:2], density = TRUE)
    upper.bounds <-  apply(rbind(apply(density.i$x, 2, max), apply(density.j$x, 2, max)), 2, max)
    lower.bounds <- apply(rbind(apply(density.i$x, 2, min), apply(density.j$x, 2, min)), 2, min)
    XYClust.Dist.Hel[as.numeric(i), as.numeric(j)] <- hellingerDistance(density.i, density.j,
                                                                        lower = lower.bounds,
                                                                        upper = upper.bounds,
                                                                        tol = 1e-6)
  }
}
XYClust.Dist.Hel.Norm <- abs(XYClust.Dist.Hel) / max(abs(XYClust.Dist.Hel))

XYClust.Dist.SKL <- array(0, dim = c(length(X.clust.labels), length(Y.clust.labels)))
for(i in X.clust.labels){
  i.indexes <- which(data.X$id == i)
  density.i <- ks::kde(data.X[i.indexes, 1:2], density = TRUE) 
  for(j in Y.clust.labels){
    j.indexes <- which(data.Y$id == j)
    density.j <- ks::kde(data.Y[j.indexes, 1:2], density = TRUE)
    upper.bounds <-  apply(rbind(apply(density.i$x, 2, max), apply(density.j$x, 2, max)), 2, max)
    lower.bounds <- apply(rbind(apply(density.i$x, 2, min), apply(density.j$x, 2, min)), 2, min)
    XYClust.Dist.SKL[as.numeric(i), as.numeric(j)] <- symmetricKLDivergence(density.i, density.j,
                                                                        lower = lower.bounds,
                                                                        upper = upper.bounds,
                                                                        tol = 1e-7, tol0 = 1e-7)
  }
}
XYClust.Dist.SKL.Norm <- XYClust.Dist.SKL / XYClust.Dist.SKL[1,1] 

par(mfrow = c(1, 1), mai = c(1.02, 0.92, 0.82, 0.42))

plot(1:4, XYClust.Dist.Norm[, 1], pch = "1", col = 1, ylim = c(0, 1),
     xlab = "Origin Cluster (Cytometry)", ylab = "Normalised Distance", lab = c(3, 5 ,7),
     main = "Distance comparison between clusters", cex = 2,
     cex.axis = 1.5, cex.lab = 1.7, cex.main = 2)
points(1:4, XYClust.Dist.Norm[, 2], pch = "2", col = 1, cex = 2)
points(1:4, XYClust.Dist.Norm[, 3], pch = "3", col = 1, cex = 2)
points(1:4, XYClust.Dist.Norm[, 4], pch = "4", col = 1, cex = 2)

points(1:4, XYClust.Dist.MMD.vanilla.Norm[, 1], pch = "1", col = 2, cex = 2)
points(1:4, XYClust.Dist.MMD.vanilla.Norm[, 2], pch = "2", col = 2, cex = 2)
points(1:4, XYClust.Dist.MMD.vanilla.Norm[, 3], pch = "3", col = 2, cex = 2)
points(1:4, XYClust.Dist.MMD.vanilla.Norm[, 4], pch = "4", col = 2, cex = 2)

points(1:4, XYClust.Dist.MMD.gauss.Norm[, 1], pch = "1", col = 3, cex = 2)
points(1:4, XYClust.Dist.MMD.gauss.Norm[, 2], pch = "2", col = 3, cex = 2)
points(1:4, XYClust.Dist.MMD.gauss.Norm[, 3], pch = "3", col = 3, cex = 2)
points(1:4, XYClust.Dist.MMD.gauss.Norm[, 4], pch = "4", col = 3, cex = 2)

points(1:4, XYClust.Dist.Hel.Norm[, 1], pch = "1", col = 4, cex = 2)
points(1:4, XYClust.Dist.Hel.Norm[, 2], pch = "2", col = 4, cex = 2)
points(1:4, XYClust.Dist.Hel.Norm[, 3], pch = "3", col = 4, cex = 2)
points(1:4, XYClust.Dist.Hel.Norm[, 4], pch = "4", col = 4, cex = 2)

points(1:4, XYClust.Dist.SKL.Norm[, 1], pch = "1", col = 5, cex = 2)
points(1:4, XYClust.Dist.SKL.Norm[, 2], pch = "2", col = 5, cex = 2)
points(1:4, XYClust.Dist.SKL.Norm[, 3], pch = "3", col = 5, cex = 2)
points(1:4, XYClust.Dist.SKL.Norm[, 4], pch = "4", col = 5, cex = 2)

points(1:4, XYClust.Dist.FR.Norm[, 1], pch = "1", col = 6, cex = 2)
points(1:4, XYClust.Dist.FR.Norm[, 2], pch = "2", col = 6, cex = 2)
points(1:4, XYClust.Dist.FR.Norm[, 3], pch = "3", col = 6, cex = 2)
points(1:4, XYClust.Dist.FR.Norm[, 4], pch = "4", col = 6, cex = 2)

plot.new()
legend("center", legend = c("", "", "", ""),
       pch = as.character(1:4), title = "Destination cluster (Cytometry 2)", horiz = TRUE, cex = 1.3)

legend("bottom", legend = c("2-Wasserstein", "MMD Vanilla", "MMD Gaussian", "Hellinger", "SKL", "FR"),
       fill = 1:6, title = "Discrepancy Measure", horiz = TRUE, cex = 1.3)
