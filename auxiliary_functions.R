hellingerDistance <- function(density1, density2, lowerLimit, upperLimit, tol) {
  return(1 - cubature::hcubature(functionToIntegrateHell, lowerLimit, upperLimit,
                       densityEstimation1 = density1, densityEstimation2 = density2,
                       tol = tol, vectorInterface = TRUE)$integral)  
}
symmetricKLDivergence <- function(density1, density2, lowerLimit, upperLimit, tol, tol0 = 1e-7) {
  return(cubature::hcubature(functionToIntegrateKL, lowerLimit, upperLimit,
                       densityEstimation1 = density1, densityEstimation2 = density2,
                       tol = tol, tol0 = tol0, vectorInterface = TRUE)$integral)  
}

functionToIntegrateHell <- function(X, densityEstimation1, densityEstimation2) {
  X <- t(X)
  matrix(sqrt(abs(predict(densityEstimation1, x = X) * predict(densityEstimation2,
                                                               x = X))),
         ncol = nrow(X))
}
functionToIntegrateKL <- function(X, densityEstimation1, densityEstimation2, tol0 = 10^(-6)) {
  X <- t(X)
  f.density.1 <- abs(predict(densityEstimation1, x = X))
  f.density.2 <- abs(predict(densityEstimation2, x = X))
  index.zeros.1 <- which(f.density.1 == 0, arr.ind = TRUE)
  index.zeros.2 <- which(f.density.2 == 0, arr.ind = TRUE)
  n1 <- length(index.zeros.1)
  n2 <- length(index.zeros.2)
  
  if(n1 > 0){
    if(n2 > 0){
      if(sum(f.density.2[index.zeros.1] < tol0) / n1 == 1 & sum(f.density.1[index.zeros.2] < tol0) / n2 == 1){
        SKLMatrix <- matrix(f.density.1 * log(f.density.1 / f.density.2) + 
                              f.density.2 * log(f.density.2 / f.density.1),
                            ncol = nrow(X))
        indexes.nan <- which(is.nan(SKLMatrix), arr.ind = TRUE)
        number.nans <- dim(indexes.nan)[1]
        SKLMatrix[indexes.nan] <- rep(0, number.nans)
      } else {
        SKLMatrix <- matrix(f.density.1 * log(f.density.1 / f.density.2) + 
                              f.density.2 * log(f.density.2 / f.density.1),
                            ncol = nrow(X))
      }
    } else {
      if(sum(f.density.2[index.zeros.1] < tol0) / n1 == 1){
        SKLMatrix <- matrix(f.density.1 * log(f.density.1 / f.density.2) + 
                              f.density.2 * log(f.density.2 / f.density.1),
                            ncol = nrow(X))
        indexes.nan <- which(is.nan(SKLMatrix), arr.ind = TRUE)
        number.nans <- dim(indexes.nan)[1]
        SKLMatrix[indexes.nan] <- rep(0, number.nans)
      } else {
        SKLMatrix <- matrix(f.density.1 * log(f.density.1 / f.density.2) + 
                              f.density.2 * log(f.density.2 / f.density.1),
                            ncol = nrow(X))
      }
    }
  } else {
    if(n2 > 0){
      if(sum(f.density.1[index.zeros.2] < tol0) / n2 == 1){
        SKLMatrix <- matrix(f.density.1 * log(f.density.1 / f.density.2) + 
                              f.density.2 * log(f.density.2 / f.density.1),
                            ncol = nrow(X))
        indexes.nan <- which(is.nan(SKLMatrix), arr.ind = TRUE)
        number.nans <- dim(indexes.nan)[1]
        SKLMatrix[indexes.nan] <- rep(0, number.nans)
      } else {
        SKLMatrix <- matrix(f.density.1 * log(f.density.1 / f.density.2) + 
                              f.density.2 * log(f.density.2 / f.density.1),
                            ncol = nrow(X))
      }
    } else {
      SKLMatrix <- matrix(f.density.1 * log(f.density.1 / f.density.2) + 
                            f.density.2 * log(f.density.2 / f.density.1),
                          ncol = nrow(X))
    }
  }
  return(SKLMatrix)
}
