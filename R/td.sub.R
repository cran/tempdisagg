SubRegressionBased <- function(y_l, X,
  conversion="sum", method="chow-lin-minrss", fr=4, vcov = "ecotrim", 
  neg.rho = TRUE, tol=1e-16, lower=0) {
  # performs temporal disaggregation for regression based methods
  #
  # Args:
  #   y_l:          vector of the low-frequency left-hand side variable
  #   X:            matrix of high-frequency indicators
  #   conversion:   type of conversion ("sum", "average", "first", "last")
  #   method:       method
  #   fr:           ratio of high-frequency units per low-frequency unit
  #   tol:          desired accuracy, passed on to optim()
  #   lower:        scalar indicating the lower limit of the parameter
  #
  # Returns:
  #   a list containing the following variables
  #     vcov            type of the variance-covariance approx.
  #     fitted.values   interpolated (and extrapolated) high frequency series
  #     p               preliminary high frequency series
  #     residuals       low-frequency residuals
  #     rho             autoregressive parameter
  #     coefficients    a named vector of coefficients
  #     se              standard errors of the coefficients
  #     s_2             ML-estimator of the variance of the high-freq. res.
  #     s_2_gls         GLS-estimator of the variance of the high-freq. res.
  #     tss             total sum of square
  #     rss             (low frequency) residual sum of square, weighted by the
  #                       variance-covariance matrix
  #     logl            log-likelihood
  #     rank            number of right hand variables (including intercept)
  #     df              degrees of freedom
  

  z <- list()
  
  # dimensions of y_l and X
  n_l <- length(y_l)
  n <- dim(X)[1]
  m <- dim(X)[2]

  # conversion matrix expanded with zeros
  C <- CalcC(n_l, conversion, fr, n)

  pm <- CalcPowerMatrix(n)

  # method specific variance-covariance functions
  if (method=="chow-lin-maxlog"){
    Objective <- function(rho){-CalcLogL(CalcQ(rho, pm), C=C, y_l=y_l, X=X)}
  } else if (method=="chow-lin-minrss"){
      z$vcov <- vcov  # append to output
      if (vcov == "ecotrim"){
        Objective <- function(rho){CalcRSS(CalcR(rho, pm), C=C, y_l=y_l, X=X)}
      } else if (vcov == "quilis"){
        Objective <- function(rho){CalcRSS(CalcQ(rho, pm), C=C, y_l=y_l, X=X)}
      } else {
        stop("wrong optim specification")
      }
  } else if (method=="litterman-maxlog"){
    Objective <- function(rho){-CalcLogL(CalcQ_Lit(X, rho), C=C, y_l=y_l, X=X)}
  } else if (method=="litterman-minrss"){      
    if (vcov == "ecotrim"){
      Objective <- function(rho){CalcRSS(CalcR(rho, pm), C=C, y_l=y_l, X=X)}
    } else if (vcov == "quilis"){
      Objective <- function(rho){CalcRSS(CalcQ(rho, pm), C=C, y_l=y_l, X=X)}
    } else {
      stop("wrong optim specification")
    }
 }

  # finding the optimal rho parameter
  if (method %in% c("chow-lin-maxlog", "chow-lin-minrss", 
      "litterman-maxlog", "litterman-minrss")){
    optimize.results <- optimize(Objective, lower=lower, upper=0.999, tol=tol,
                            maximum=FALSE)
    rho <- optimize.results$minimum
    if (rho < 0){
      if (neg.rho){
        warning("negative rho value", immediate.=TRUE)
      } else {
        rho <- 0
      }
    }
    if (rho < -0.998 | rho > 0.998){
      warning("boundary solution for rho", immediate.=TRUE)
    }
  } else if (method=="fernandez"){
    rho <- 1
  } else if (method=="ols"){
    rho <- 0
  }

  Q       <- CalcQ(rho, pm)
  b       <- CalcBeta(Q, C, y_l, X)
  s_2     <- CalcSigma2(Q, C, y_l, X)
  s_2_gls <- s_2 * n_l / (n_l - m) 
  S       <- s_2 * Q
  CQC_inv <- solve(C %*% Q %*% t(C))
  CSC_inv <- CQC_inv /s_2
  p       <- as.numeric(X %*% b)
  se      <- as.numeric(sqrt(diag(s_2_gls * solve(t(X) %*% t(C) %*% CQC_inv 
                                                  %*% C %*% X))))

  # total sum of squares
  e <- matrix(data=1, nrow=n_l, ncol=1)
  y_l_bar <- as.numeric(t(e) %*% CSC_inv %*% y_l / t(e) %*% CSC_inv %*% e)
  tss  <- as.numeric(t(y_l - y_l_bar) %*% CSC_inv %*% (y_l - y_l_bar))
  rss  <- CalcRSS(S, C, y_l, X)
  logl <- CalcLogL(Q, C, y_l, X)

  # distribution matrix
  D <- Q %*% t(C) %*% solve(C %*% Q %*% t(C))

  # low frequency residuals
  u_l <- as.numeric(y_l - C %*% p)

  # final series
  y <- as.numeric(p + D %*% u_l)
 
  # output
  z$fitted.values    <- y
  z$p                <- p
  z$residuals        <- u_l
  z$rho              <- rho
  z$coefficients     <- b
  z$se               <- se
  z$s_2              <- s_2
  z$s_2_gls          <- s_2_gls
  z$tss              <- tss
  z$rss              <- rss
  z$logl             <- logl
  z$rank             <- m
  z$df               <- n_l - m   
  z
}


SubDenton <- function(y_l, X, conversion, method, fr, 
                      criterion="proportional", h=1) {
  # performs temporal disaggregation for denton methods
  #
  # Args:
  #   y_l:          vector of the low-frequency left-hand side variable
  #   X:            matrix of high-frequency indicators
  #   conversion:   type of conversion ("sum", "average", "first", "last")
  #   method:       method
  #   fr:           ratio of high-frequency units per low-frequency unit
  #
  # Returns:
  #   a list containing the following variables
  #     fitted.values   interpolated (and extrapolated) high frequency series
  #     p               preliminary high frequency series
  #     residuals       low-frequency residuals
  
  if (dim(as.matrix(X))[2] > 1){
    stop("Right hand side is not a vector, only one series allowed in
         Denton methods")
  }
  if (!(criterion %in% c("additive", "proportional"))) {
    stop("criterion for Denton methods must be additive or proportional")
  }
  
  # uniform is a special case of denton
  if (method == "uniform"){
    h <- 0
    criterion <- "additive"
    method <- "denton"
  }
  
  # dimensions of y_l and X
  n_l <- length(y_l)
  n <- length(as.numeric(X))

  # conversion matrix expanded with zeros
  C <- CalcC(n_l, conversion, fr, n)

  D <- diag(n)
  diag(D[2:n, 1:(n-1)]) <- -1
  D_0 <- diag(n)
  X_inv <- diag(1 / (as.numeric(X)/mean(X)))

  if (h==0) {
    if(criterion=="proportional") {
      D_0 <- D_0 %*% X_inv
    }
    D_1 <- D_0
  } else if (h>0) {
    for (i in 1:h) {
      D_0 <- D%*%D_0
    }
    if(criterion=="proportional") {
      D_0 <- D_0 %*% X_inv
    }
  } else stop("wrong specification of h")

  # low frequency residuals
  u_l <- as.numeric(y_l - C %*% X)
  
  if (method=="denton-cholette"){
    
    D_1 <- D_0[-(1:h),]
    A <- t(D_1)%*% D_1
    
    # Eq. (2.2) from Denton (1971); Eq (6.8) from Cholette and Dagum (2006)
    y <- solve(
      rbind(cbind(A, t(C)), cbind(C, matrix(0, nrow=n_l, ncol=n_l)))
      ) %*% rbind(
      cbind(A, matrix(0, nrow=n, ncol=n_l)),
      cbind(C, diag(1, nrow=n_l, ncol=n_l))
      ) %*% matrix(c(X, u_l))
    
    # final series
    y <- y[1:n]
    
  } else if (method=="denton"){
    
    D_1 <- D_0
    
    # Denton (1971), in the text below Eq. (2.2)
    Q <- solve(t(D_1) %*% D_1)
    
    # distribution matrix
    D <- Q %*% t(C) %*% solve(C %*% Q %*% t(C))
  
    # final series
    y <- as.numeric(X + D %*% u_l)
  }
  
  # output
  z <- list()
  z$fitted.values <- y
  z$p             <- X
  z$residuals     <- u_l
  z$criterion     <- criterion
  z$h             <- h

  z
}


