require(optimalFlowData)
require(reshape2)
require(ggplot2)
require(tclust)
require(flowMap)
require(ks)

database <- optimalFlowData::buildDatabase(
  dataset_names = paste0('Cytometry', c(2:5, 7:9, 12:17, 19, 21)),
  population_ids = c('Monocytes', 'CD4+CD8-', 'Mature SIg Kappa', 'TCRgd-'))
pairs(database[[1]][, c(4, 3, 9)], col = droplevels(database[[1]][, 11]))
plot(database[[1]][, c(4, 3)], col = droplevels(database[[1]][, 11]), 
     xlim = c(0, 8000), ylim = c(0,8000))
points(database[[2]][, c(4, 3)], col = droplevels(database[[2]][, 11]), pch = 2)
points(database[[3]][, c(4, 3)], col = droplevels(database[[3]][, 11]), pch = 3)
points(database[[4]][, c(4, 3)], col = droplevels(database[[4]][, 11]), pch = 4)
points(database[[5]][, c(4, 3)], col = droplevels(database[[5]][, 11]), pch = 4)
plot(database[[3]][, c(4, 3)], col = droplevels(database[[3]][, 11]), pch = 3)

X <- database[[1]][, c(4, 3)]
X.labels <- database[[1]][, 11]
Y <- database[[5]][, c(4, 3)]
Y.labels <- database[[5]][, 11]

plot(X, col = "blue", main = "Two ungated cytometry datasets")
points(Y, col = "orange", pch = 4)
legend("bottomleft", legend = c("Cytometry 1", "Cytometry 2"), pch = c(1,4), col = c("blue", "orange"), bty = "n")

write.csv(X, file = "data1.csv", row.names = FALSE)
write.csv(Y, file = "data2.csv", row.names = FALSE)

density.Y <- ks::kde(Y, density = TRUE) 
density.X <- ks::kde(X, density = TRUE) 

plot(density.Y, display = "filled.contour", cont = seq(5, 95, 5),
     xlim = c(0, 7000), ylim = c(0, 8000), alpha = 0.55,
     drawlabels = TRUE,
     main = "Density estimate of two ungated cytometries")
plot(density.X, display = "filled.contour2", cont = seq(5, 95, 5),
     xlim = c(0, 7000), ylim = c(0, 8000), add = TRUE, alpha = 0.55)

t0 <- Sys.time()
FRdistXY <- flowMap::getFR(X, Y)
t1 <- Sys.time()
print(t1 - t0)

X.clustered <- tclust::tclust(X, k = 4, alpha = 0)
Y.clustered <- tclust::tclust(Y, k = 4, alpha = 0)

plot(X, col = "blue", pch = as.character(X.clustered$cluster), main = "Two automatically-gated cytometry datasets")
points(Y, col = "orange", pch = as.character(Y.clustered$cluster))
legend("bottomleft", legend = c("Cytometry 1", "Cytometry 2"),
       fill = c("blue", "orange"), bty = "n")

pooled.barycenter <- rbind(X, Y)
plot(pooled.barycenter)

N.barr <- 1.5e+3
w1 <- 0.5
w2 <- 0.5
n.X <- dim(X)[1]
n.Y <- dim(Y)[1]
X.bar.index <- sample(1:n.X, round(w1 * N.barr))
Y.bar.index <- sample(1:n.Y, round(w2 * N.barr))

XY.MMD.bar <- rbind(X[X.bar.index,], Y[Y.bar.index, ])

XY.Pool.bar <- rbind(X, Y)[sample(1:(n.X + n.Y), N.barr), ]

XY.Wasser.Bar <- read.csv("~/pCloudDrive/R/Capitulo_Libro_Cytometry/wasser_bar.csv", header=FALSE)
colnames(XY.Wasser.Bar) <- colnames(X)

par(mfrow = c(1, 1), mai = c(1.02, 0.92, 0.82, 0.42))
plot(XY.Wasser.Bar, xlim = c(0, 7000), ylim = c(0, 8000), main = "Ungated datasets: barycenters with 1500 points",
     cex = 1.5,
     cex.main = 2, cex.axis = 1.5, cex.lab = 1.7)
points(XY.MMD.bar, col = 2, pch = 2)
points(XY.Pool.bar, col = 3, pch = 3)
legend("bottomleft", legend = c("Wasserstein", "MMD", "Pooling"),
       pch = 1:3, col = 1:3, bty = "n", cex = 1.5)

avrg.prop <- colMeans(rbind(table(X.clustered$cluster) / sum(table(X.clustered$cluster)),
                            table(Y.clustered$cluster) / sum(table(Y.clustered$cluster))))
avrg.prop <- avrg.prop / sum(avrg.prop)

clust.bar.size <- round(avrg.prop * N.barr)

write.csv(X[sample(which(X.clustered$cluster == 1), 4000), ], file = "data1C1.csv", row.names = FALSE)
write.csv(Y[sample(which(Y.clustered$cluster == 1), 1000), ], file = "data2C1.csv", row.names = FALSE)

write.csv(X[sample(which(X.clustered$cluster == 2), 4000), ], file = "data1C2.csv", row.names = FALSE)
write.csv(Y[sample(which(Y.clustered$cluster == 2), 1000), ], file = "data2C2.csv", row.names = FALSE)

write.csv(X[X.clustered$cluster == 3, ], file = "data1C3.csv", row.names = FALSE)
write.csv(Y[Y.clustered$cluster == 3, ], file = "data2C3.csv", row.names = FALSE)

write.csv(X[X.clustered$cluster == 4, ], file = "data1C4.csv", row.names = FALSE)
write.csv(Y[Y.clustered$cluster == 4, ], file = "data2C4.csv", row.names = FALSE)

XY.Wasser.Bar.2.1 <- read.csv("~/pCloudDrive/R/Capitulo_Libro_Cytometry/wasser_bar_C1.csv", header=FALSE)
XY.Wasser.Bar.2.2 <- read.csv("~/pCloudDrive/R/Capitulo_Libro_Cytometry/wasser_bar_C2.csv", header=FALSE)
XY.Wasser.Bar.2.3 <- read.csv("~/pCloudDrive/R/Capitulo_Libro_Cytometry/wasser_bar_C3.csv", header=FALSE)
XY.Wasser.Bar.2.4 <- read.csv("~/pCloudDrive/R/Capitulo_Libro_Cytometry/wasser_bar_C4.csv", header=FALSE)
XY.Wasser.Bar.2 <- rbind(XY.Wasser.Bar.2.1, XY.Wasser.Bar.2.2,
                         XY.Wasser.Bar.2.3, XY.Wasser.Bar.2.4)
colnames(XY.Wasser.Bar.2) <- colnames(X)


XY.Pool.bar.2.1 <- rbind(X[sample(which(X.clustered$cluster == 1), round(clust.bar.size[1]) * (1 - n.Y/n.X)), ],
                         Y[sample(which(Y.clustered$cluster == 1), round(clust.bar.size[1]) * (n.Y/n.X)), ])
XY.Pool.bar.2.2 <- rbind(X[sample(which(X.clustered$cluster == 2), round(clust.bar.size[2]) * (1 - n.Y/n.X)), ],
                         Y[sample(which(Y.clustered$cluster == 2), round(clust.bar.size[2]) * (n.Y/n.X)), ])
XY.Pool.bar.2.3 <- rbind(X[sample(which(X.clustered$cluster == 3), round(clust.bar.size[3]) * (1 - n.Y/n.X)), ],
                         Y[sample(which(Y.clustered$cluster == 3), round(clust.bar.size[3]) * (n.Y/n.X)), ])
XY.Pool.bar.2.4 <- rbind(X[sample(which(X.clustered$cluster == 4), round(clust.bar.size[4]) * (1 - n.Y/n.X)), ],
                         Y[sample(which(Y.clustered$cluster == 4), round(clust.bar.size[4]) * (n.Y/n.X)), ])
XY.Pool.bar.2 <- rbind(XY.Pool.bar.2.1, XY.Pool.bar.2.2, 
                       XY.Pool.bar.2.3, XY.Pool.bar.2.4)

XY.MMD.bar.2.1 <- rbind(X[sample(which(X.clustered$cluster == 1), round(clust.bar.size[1]) * 0.5), ],
                         Y[sample(which(Y.clustered$cluster == 1), round(clust.bar.size[1]) * 0.5), ])
XY.MMD.bar.2.2 <- rbind(X[sample(which(X.clustered$cluster == 2), round(clust.bar.size[2]) * (1 - 0.5)), ],
                         Y[sample(which(Y.clustered$cluster == 2), round(clust.bar.size[2]) * (0.5)), ])
XY.MMD.bar.2.3 <- rbind(X[sample(which(X.clustered$cluster == 3), round(clust.bar.size[3]) * (1 - 0.5)), ],
                         Y[sample(which(Y.clustered$cluster == 3), round(clust.bar.size[3]) * (0.5)), ])
XY.MMD.bar.2.4 <- rbind(X[sample(which(X.clustered$cluster == 4), round(clust.bar.size[4]) * (1 - 0.5)), ],
                         Y[sample(which(Y.clustered$cluster == 4), round(clust.bar.size[4]) * (0.5)), ])
XY.MMD.bar.2 <- rbind(XY.MMD.bar.2.1, XY.MMD.bar.2.2, 
                       XY.MMD.bar.2.3, XY.MMD.bar.2.4)

plot(XY.Wasser.Bar.2, col = 1,
     pch = c(rep("1", dim(XY.Wasser.Bar.2.1)[1]),
             rep("2", dim(XY.Wasser.Bar.2.2)[1]),
             rep("3", dim(XY.Wasser.Bar.2.3)[1]),
             rep("4", dim(XY.Wasser.Bar.2.4)[1])),
     xlim = c(0, 7000), ylim = c(0, 8000), main = "Gated datasets: barycenters with 1500 points",
     cex = 1.5,
     cex.main = 2, cex.axis = 1.5, cex.lab = 1.7)
points(XY.MMD.bar.2, col = 2,
       pch = c(rep("1", dim(XY.MMD.bar.2.1)[1]),
               rep("2", dim(XY.MMD.bar.2.2)[1]),
               rep("3", dim(XY.MMD.bar.2.3)[1]),
               rep("4", dim(XY.MMD.bar.2.4)[1])),
       cex = 1.5)
points(XY.Pool.bar.2, col = 3,
       pch = c(rep("1", dim(XY.Pool.bar.2.1)[1]),
               rep("2", dim(XY.Pool.bar.2.2)[1]),
               rep("3", dim(XY.Pool.bar.2.3)[1]),
               rep("4", dim(XY.Pool.bar.2.4)[1])),
       cex = 1.5)
legend("bottomleft", legend = c("Wasserstein", "MMD", "Pooling"),
       pch = 1:3, col = 1:3, bty = "n", cex = 1.5)

plot(XY.Wasser.Bar.2, col = 1,
     xlim = c(0, 7000), ylim = c(0, 8000), main = "Gated datasets: barycenters with 1500 points",
     cex = 1.5,
     cex.main = 2, cex.axis = 1.5, cex.lab = 1.7)
points(XY.MMD.bar.2, col = 2,
       pch = 2)
points(XY.Pool.bar.2, col = 3,
       pch = 3)
legend("bottomleft", legend = c("Wasserstein", "MMD", "Pooling"),
       pch = 1:3, col = 1:3, bty = "n", cex = 1.5)