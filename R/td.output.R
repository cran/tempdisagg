print.td <- function(x, ...){
  # calls print.lm
  #
  # Args:
  #   x:           an object of class "td"
  #
  # Returns: 
  #   prints the object with print.lm as a side effect
  print.lm(x, ...)
  
  # for the internal version only:
  cat("\nUse summary() for details. \nUse predict() to extract the final series.
      \nUse ?td to see the help file.\n")
}

summary.td <- function(object, ...){
  # prepares a summary for an object of class "td"
  #
  # Args:
  #   object: an object of class "td"
  #
  # Returns: 
  #   an object of class "summary.td"

  # change name from 'object' (default of the generic function) to standard 'x'
  x <- object
  
  # build output on top of the input
  z <- object

  # number of observations and number of indicator series
  z$n_l <- length(x$actual)
  if (is.null(dim(x$model))){
    m <- 1
    z$n <- length(x$model)
  } else {
    m <- dim(x$model)[2]
    z$n <- dim(x$model)[1]
  }

  # derivative statististics
  if (!is.null(x$rss)){
    z$sigma          <- sqrt(x$rss / z$df)
    z$r.squared      <- 1 - x$rss / x$tss
    z$adj.r.squared  <- 1 - (x$rss * (z$n_l-1)) / (x$tss * (z$n_l - m))
    z$aic            <- log(x$rss / z$n_l) + 2 * (m/z$n_l)
    z$bic            <- log(x$rss / z$n_l) + log(z$n_l) * (m / z$n_l)
  }
  z$ar_l           <- cor(z$residuals[-1], z$residuals[-length(z$residuals)])

  # coefficents matrix
  if (!is.null(coef(x))){
    est  <- coef(x)
    se   <- x$se
    tval <- est/se
    pval <- 2 * pnorm(-abs(tval))
    
    z$coefficients <- cbind(est, se, tval, 
                            2 * pt(abs(tval), z$df, lower.tail = FALSE))
    dimnames(z$coefficients) <- list(names(est),
                                     c("Estimate", "Std. Error", "t value", 
                                       "Pr(>|t|)"))
  }
  class(z) <- "summary.td"
  z
}

print.summary.td <- function (x, digits = max(3, getOption("digits") - 3), 
      signif.stars = getOption("show.signif.stars"), ...) {
  # prints a summary for an object of class "td"
  #
  # Args:
  #   x:           an object of class "summary.td"
  #   digits:      the number of significant digits to use when printing.
  #   signif.stars whether to show significance stars in output
  #
  # Returns: 
  #   prints the summary as its side effect
  
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "")
  resid <- x$residuals
  cat("Residuals:\n")
  if (length(resid) > 5) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    quantile.resid <- zapsmall(quantile(resid), digits + 1)
    print(structure(quantile.resid, names = nam), digits = digits)
  } else {
    print(resid, digits = digits)
  }
  if (is.null(coef(x))) {
    cat("\nNo Coefficients\n")
  } else {
      cat("\nCoefficients:\n")
      coefs <- coef(x)
      if (!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4, 
                        dimnames = list(cn, colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
      }
      printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
        na.print = "NA")
  }
  
  cat("\n'", x$method, "' disaggregation with '", x$conversion, 
      "' conversion\n", sep = "")
  cat(x$n_l, "low-freq. obs. converted to", x$n, "high-freq. obs.")
  if (!is.null(x$adj.r.squared)) {
    cat("\nAdjusted R-squared:", formatC(x$adj.r.squared, digits = digits))
  }
  if (!is.null(x$rho)) {
    cat("\tAR1-Parameter:", formatC(x$rho, digits = digits))
  }
  if (!is.null(x$criterion)) {
    cat("\ncriterion:", x$criterion, "\torder of differencing 'h':", x$h)
  }
  cat("\n")
  invisible(x)
}

plot.td <- function(x, ...){
  # plot the predicted and acutal low frequency series, and residuals
  #
  # Args:
  #   object:      an object of class "td"
  #   ...          further arguments, not used
  #
  # Returns: 
  #   a two panel plot as its side effect
  
  old.par <- par(no.readonly=TRUE)  # backup par settings 
  par(mfrow=c(2,1))                    
  ts.plot(ts.intersect(x$actual, x$actual - x$residuals), 
      lty=c("solid", "dashed"), col=c("black", "red"), 
      ylab="actual vs. predicted (red)"); grid()
  ts.plot(ts.intersect(x$residuals, 0), lty=c("solid", "dashed"),
      col=c("black", "red"), ylab="residuals"); grid()
  on.exit(par(old.par))  # restore par settings
}

predict.td <- function(object, ...) {
  # extracts the final series
  #
  # Args:
  #   object:      an object of class "td"
  #   ...          further arguments, not used
  #
  # Returns: 
  #   an object of class "ts", containing the final series
  object$fitted.values
}
