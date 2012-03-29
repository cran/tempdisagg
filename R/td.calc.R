CalcC <- function(n_l, conversion, fr, n=NULL){
  # calculates the conversion matrix C, optionally expanded with zeros
  #
  # Args:
  #   n_l:          number of low-frequency observations
  #   conversion:   char string indicating the type of conversion
  #                 ("sum", "average", "first", "last")
  #   fr:           ratio of high-frequency units per low-frequency unit
  #   n:            number of high-frequency observations 
  #                 (matrix will be expanded by 0 if n is provided and > n_l*fr)
  #
  # Returns: 
  #   The conversion matrix C

  # set conversion.weights according to type of conversion
  if (conversion=="sum") {
    conversion.weights <- rep(1, fr)
  } else if (conversion=="average") {
    conversion.weights <- rep(1, fr)/fr
  } else if (conversion=="first") {
    conversion.weights <- numeric(fr)
    conversion.weights[1] <- 1
  } else if (conversion=="last") {
    conversion.weights <- numeric(fr)
    conversion.weights[fr] <- 1  
  } else stop("Wrong type of conversion")

  # compute the conversion matrix
  C <- kronecker(diag(n_l), t(conversion.weights))
  if(!is.null(n)){
    C <- cbind(C, matrix(0, nrow=n_l, ncol=n - fr * n_l))
  }
  C
}

CalcBeta <- function(Q, C, y_l, X){
  # calculates GLS/ML-estimator of beta
  #
  # Args:
  #   Q:           s_2-factored-out VCov (S = s_2 * Q)
  #                (instead of Q, S and R can be used to estimate b)
  #   S:           variance-covariance matrix (VCov)
  #   R:           correlation matrix
  #                b is not affected by whether S, Q or R is used
  #   C:           conversion matrix, with or without zeros
  #   y_l:         vector of the low-frequency left-hand side variable
  #   X:           matrix of high-frequency indicators
  #
  # Returns: 
  #   b:           GLS/ML-estimator of beta, a vector

  CQC_inv <- solve(C %*% Q %*% t(C))
  as.numeric(solve(t(X) %*% t(C) %*% CQC_inv %*% C %*% X) %*% t(X) %*% t(C) 
             %*% CQC_inv %*% y_l)
}

CalcSigma2 <- function(Q, C, y_l, X){
  # calculates s_2, the ML-estimator of the variance of the high-freq. res.
  #
  # Args:
  #   Q:           s_2-factored-out variance-covariance matrix
  #   C:           conversion matrix, with or without zeros
  #   y_l:         vector of the low-frequency left-hand side variable
  #   X:           matrix of high-frequency indicators
  #
  # Returns: 
  #   s_2:         ML-estimator of the variance of the high-freq. res., a scalar

  n_l <- length(y_l)

  # ML-Estimator for b (equal to GLS)
  b <- CalcBeta(Q, C, y_l, X)

  # ML-Estimator for s_2
  u_l <- y_l - C %*% X %*% b
  CQC_inv <- solve(C %*% Q %*% t(C))
  as.numeric((t(u_l) %*% CQC_inv %*% u_l) / n_l)
}

CalcLogL <- function(Q, C, y_l, X){
  # calculates the (low frequency) log likelihood (LogL)
  #
  # Args:
  #   Q:           variance-covariance matrix, correlation matrix can be used
  #                in some cases.
  #   C:           conversion matrix, with or without zeros
  #   y_l:         vector of the low-frequency left-hand side variable
  #   X:           matrix of high-frequency indicators
  #
  # Returns: 
  #   logl:        log likelihood, a scalar

  n_l <- length(y_l)
  m <- dim(X)[2]

  # ML-Estimator for b (equal to GLS)
  b <- CalcBeta(Q, C, y_l, X)

  # ML-Estimator for s_2
  s_2 <- CalcSigma2(Q, C, y_l, X)

  # Log-Likelihood
  u_l <- y_l - C %*% X %*% b
  as.numeric(- n_l / 2 - n_l * log(2 * pi) / 2 - n_l * log(s_2) / 2 - 
              log(det(C %*% Q %*% t(C))) / 2)
}

CalcRSS <- function(S, C, y_l, X){
  # calculates the (low frequency) residual sum of square (RSS), 
  #   weighted by the variance-covariance matrix
  #
  # Args:
  #   S:           variance-covariance matrix, correlation matrix can be used
  #                in some cases.
  #   C:           conversion matrix, with or without zeros
  #   y_l:         vector of the low-frequency left-hand side variable
  #   X:           matrix of high-frequency indicators
  #
  # Returns: 
  #   rss:         weighted residual sum of squares, a scalar

  # ML-Estimator for b (equal to GLS)
  # (b is not affected by whether S, Q or R is used)
  b <- CalcBeta(S, C, y_l, X)
  
  u_l <- y_l - C %*% X %*% b
  CSC_inv <- solve(C %*% S %*% t(C))
  as.numeric(t(u_l) %*% CSC_inv %*% u_l)
}

CalcPowerMatrix <- function(n){
  # calculates a symetric 'power' matrix with 0 on the diagonal, 
  #   1 in the subsequent diagonal, and so on.
  #
  # Args:
  #   n:           number of high-frequency observations
  #
  # Returns: 
  #   power matrix  

  mat <- diag(n)
  abs(row(mat) - col(mat))
}

# covariance-, correlation- and in-between-matrices
CalcR <- function(rho, pm){
  rho^pm
}

CalcQ <- function(rho, pm){
  (1/(1-rho^2)) * CalcR(rho, pm)
}

CalcS <- function(rho, pm, s_2){
  s_2 * CalcQ(rho, pm)
}

CalcQ_Lit <- function(X, rho=0) {
  # calculates the (pseudo) variance-covariance matrix (and stats) for 
  #   a Random Walk (with opt. AR1)
  #
  # Args:
  #   X:            matrix of high-frequency indicators
  #   rho:          if != 0, a AR1 is added to the RW (Litterman)
  #
  # Returns:
  #   Q_Lit:        pseudo variance-covariance matrix

  # dimension of y_l
  n <- dim(X)[1]

  # calclation of S
  H <- diag(n)
  D <- diag(n)
  diag(D[2:nrow(D), 1:(ncol(D)-1)]) <- -1
  diag(H[2:nrow(H), 1:(ncol(H)-1)]) <- -rho
  Q_Lit <- solve(t(D)%*%t(H)%*%H%*%D)

  # output
  Q_Lit
}
