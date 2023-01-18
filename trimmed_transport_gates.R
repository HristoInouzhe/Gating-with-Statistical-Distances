require(transport)

# First we produce, through optimal trimming, a pair of cytometris 
# based on the originals, which are more similar to each other

# X and Y are defined in plots_figures_barycenters.R
 
n.X <- dim(X)[1]
n.Y <- dim(Y)[1]

N_0 <- n.X
N_1 <- n.Y

trimm.level.minority <- 0

alpha_1 <- (1 - N_1 / N_0) + trimm.level.minority
alpha_2 <- trimm.level.minority

P <- c(rep(1 / ((1 - alpha_1) * N_0), N_0), alpha_2 / (1 - alpha_2))
Q <- c(rep(1 / ((1 - alpha_2) * N_1), N_1), alpha_1 / (1 - alpha_1))

# This step is taxing on RAM memory since it saves on memory a matrix of size n.X * n.Y
t0 <- Sys.time()
distance_matrix <- (as.matrix(dist(rbind(X, Y)))^2)[1:N_0, (N_0 + 1): (N_0 + N_1)]
t1 <- Sys.time()
print(t1 - t0)

cost.matrix <- distance_matrix

cost.matrix <- cbind(cost.matrix, rep(0, N_0))
cost.matrix <- rbind(cost.matrix,
                     c(rep(0, N_1), ((1 - alpha_1) / (alpha_1 * (1 - alpha_2))) * sum(cost.matrix)))

t0 <- Sys.time()
trimmed.transport <- transport::transport(P, Q, cost.matrix)
t1 <- Sys.time()
print(t1 - t0)

trimmed.transport_1 <- trimmed.transport[trimmed.transport$from != (N_0 + 1), ]
trimmed.transport_1 <- trimmed.transport_1[trimmed.transport_1$to != (N_1 + 1), ]

non.trimmed.0 <- unique(trimmed.transport_1$from)
non.trimmed.1 <- unique(trimmed.transport_1$to)

N_0_1 <- length(non.trimmed.0)
N_1_1 <- length(non.trimmed.1)


# We define the new cytometries as subsets of the original ones

origin <- X[non.trimmed.0, ]
destination <- Y[non.trimmed.1, ]

# We do the ploting for both markers

plot(density(origin[, 1]),
     main = "Gate transportation", xlab = "CD4+CD20:PB-A LOGICAL", lwd = 2,
     cex.axis = 1.5, cex.lab = 1.7, cex.main = 2)
lines(density(destination[, 1]), col = 2, lwd = 2)
abline(v = 5250, col = 1, lty = 2, lwd = 2)
abline(v = quantile(destination[, 1],
                    ecdf(origin[, 1])(5250)),
       col = 2, lty = 2, lwd = 2)
legend("topleft", legend = c("Origin cytometry", "Destination cytometry", "Gate at origin", "Transported gate"),
       lty = c(1, 1, 2, 2), col = c(1, 2, 1, 2), 
       cex = 1.4, bty = "n")

plot(density(origin[, 2]),
     main = "Gate transportation", xlab = "CD3:APC-A LOGICAL", lwd = 2,
     cex.axis = 1.5, cex.lab = 1.7, cex.main = 2)
lines(density(destination[, 2]), col = 2, lwd = 2)
abline(v = 5000, col = 1, lty = 2, lwd = 2)

# Next is the transport of the original gate,
# 5000, to its transported counterpart 

abline(v = quantile(destination[, 2],
                    ecdf(origin[, 2])(5000)),
       col = 2, lty = 2, lwd = 2)
legend("topleft", legend = c("Origin cytometry", "Destination cytometry", "Gate at origin", "Transported gate"),
       lty = c(1, 1, 2, 2), col = c(1, 2, 1, 2), 
       cex = 1.4, bty = "n")
