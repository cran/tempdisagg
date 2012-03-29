td <- function(formula, conversion="sum", method="chow-lin-minrss", 
               to="quarterly", start=NULL, end=NULL, ...) {
  # performs a temporal disagregation
  #
  # Args:
  #   formula:      a formula specifying the temporal disagregation model
  #   conversion:   type of conversion ("sum", "average", "first", "last")
  #   method:       method of temporal disagregation
  #   to:           destination frequency ("quarterly" or "monthly"), or scalar 
  #                 only for Denton without indicators
  #   start:        optional start point of the disagregation (final series
  #                 will be shortened)
  #   end:          optional end point
  #   ...           Further arguments, passed on to the subfunctions
  #
  # Returns:
  #   An object of class "td", containing the components the subroutines and the 
  #     the following objects:
  #   method        method of temporal disaggregation
  #   call          function call
  #   fr            the ratio between the low- and high-frequency series
  #   conversion    type of conversion
  #   actual        actual values of the low-frequeny series
  #   model         a matrix containing the indicators (and a intercept if
  #                   present)
  #
  # Remarks:
  #   This function deals with the formula interface, the time-series properties
  #   and allows for optional shortening of the series. The estimation itself is
  #   done by subfunctions stating with Sub...

  cl <- match.call()
  # check input consistency
  # dont allow length=2 vectors as start or end inputs
  if (any(c(length(start), length(end)) > 1)){
    stop("'start' or 'end' must be specified as a decimal fraction")
  }

  # prepare Formula, extract names and data
  # extract X (right hand side, high frequency) formula, names and data
  X.formula <- formula; X.formula[[2]] <- NULL
  X.series.names <- all.vars(X.formula)

  # extract y_l (left hand side, low frequency) formula, values and names
  y_l.formula <- formula[[2]]
  y_l.series <- eval(y_l.formula, envir=environment(formula))
  y_l.name <- deparse(y_l.formula)

  # set ts.mode
  # 1. is y_l.series a time series? if yes, set ts.mode to TRUE
  if (is.ts(y_l.series)){
    ts.mode <- TRUE
    # 2. is there a X? is it a time series? if not, set ts.mode to FALSE
    if (length(X.series.names) > 0) {
      if (!is.ts(get(X.series.names[1], envir=environment(X.formula)))){
        ts.mode <- FALSE
      }
    }
  } else{
    ts.mode <- FALSE
  }

  # optionally using the time series attributes
  if (ts.mode) {
    # Preparing the y_l.series
    if (is.null(start)) {
      start <- time(y_l.series)[!is.na(y_l.series)][1]
    }
    f_l <- frequency(y_l.series)

    if (length(X.series.names)>0){  # Preparing the X.series
      # extract the first series (to extract start and frequency)
      X.series.first <- eval(X.formula[[2]], envir=environment(X.formula))
      X.start <- time(X.series.first)[!is.na(X.series.first)][1]
      f <- frequency(X.series.first)
      fr <- f/f_l
      if (X.start > start){
        start <- floor(X.start) + 
          (ceiling(((X.start - floor(X.start)) * f) / fr)) / f_l
      }
      
    } else {  # If no X is specified
      if (is.numeric(to)){  # frequency specified by a number
        f <- to
      } else if (is.character(to)){  # frequency specified by a char string
        if (to=="quarterly"){
        f <- 4
        } else if (to=="monthly"){
        f <- 12
        } else {
        stop("unknown character string as the 'to' argument")
        }
      } else stop ("wrong specification of the 'to' argument")
      fr <- f/f_l
      X.start <- start
    }

    # define X.end, if 'end' is specified
    if(!is.null(end)){
      X.end <- floor(end) + (fr * ((end-floor(end)) * f_l + 1) - 1) / f
    } else {
      X.end <- NULL
    }

    # final y_l data matrix
    y_l.series <- window(y_l.series, start=start, end=end)
    y_l <- as.matrix(y_l.series)

    # final X data matrix
    # if there is one or more RHS Variables
    if (length(X.series.names) > 0){
      # final RHS data matrix, also deletes NAs
      X <- model.matrix(X.formula)

      # reduce the RHS dataset if start or end are specified
      if (ts.mode){
        X.names <- dimnames(X)[[2]]
        X <- ts(X, start=X.start, frequency=f)
        X <- window(X, start=start, end=X.end)
        # convert to matrix
        X <- matrix(X, nrow=nrow(X), ncol=ncol(X))
        dimnames(X) <- list(NULL, X.names)
      }
    } else {
      # if there is no X Variables, set it to a constant ('Denton' Methods)
      X <- rep(1, times=length(y_l.series) * fr)
      if (!(method %in% c("denton-cholette", "denton", "uniform"))) {
        stop ("No indicator specified: only works with denton,
              denton-cholette or uniform")
      }
    }
  }

  # actual estimation 
  if (method %in% c("chow-lin-maxlog", "chow-lin-minrss", 
      "litterman-maxlog", "litterman-minrss", "fernandez", "ols")){
    z <- SubRegressionBased(y_l=y_l, X=X, conversion=conversion, method=method,
                            fr=fr, ...)
  } else if (method %in% c("denton-cholette", "denton", "uniform")){
    z <- SubDenton(y_l=y_l, X=X, conversion=conversion, method=method,
                   fr=fr, ...)
  } else {
    stop("method does not exist")
  }

  # add coefficent names to output
  if (!is.null(z$coefficients)) {
    names(z$coefficients) <- X.names
    names(z$se) <- names(z$coefficients)
  }

  # additional output
  z$method             <- method
  z$call               <- cl
  z$name               <- y_l.name
  z$fr                 <- fr
  z$conversion         <- conversion
  z$actual             <- y_l.series
  z$model              <- X
  if (ts.mode) {
    z$model            <- ts(z$model, start=start, frequency=f)
    z$p                <- ts(z$p, start=start, frequency=f)
    z$fitted.values    <- ts(z$fitted.values, start=start, frequency=f)
    z$residuals        <- ts(z$residuals, start=start, frequency=f_l)
    z$actual           <- ts(z$actual, start=start, frequency=f_l)
  }
  class(z) <- "td"
  z
}


