
#' Calculate the pdf.
#'
#' @param t The
#' @param last.spike The time of the previous spike.
#' @param hyper The ISI parameter value
#' @param end.time The end time of the experiment
#' @param x The intensity function --- defined in 0,end.time
#' @param X The integral of the intensity function
#' @param ISI.type The ISI distribution
#' @param do.log Flag for whether to calculate the pdf on the log scale.
#'
#' @return The probability density of a spike at time t given the last spike was at last.spike.
#' @export
#'
#' @examples
PDF <- function(t, last.spike, hyper, end.time, x, X, ISI.type,do.log = F){
  # Check the ISI type is allowable
  ISIs <- c('Gamma', 'Exponential', 'InverseGaussian', 'LogNormal', 'Weibull',
            'InverseGaussianOLD','Gamma2', 'InverseGaussian2', 'LogNormal2', 'Weibull2' )
  if(!(ISI.type %in% ISIs)){ stop("Invalid ISI.type. See help documentation for allowable ISI.type.")}


  step.size <- end.time/(length(x)-1)
  x <- x[t/step.size]
  if(last.spike - 0 < 1e-10){
    X <- X[t/step.size]
  }
  else{
    X <- X[t/step.size] - X[last.spike/step.size]
  }
  if(!do.log){
    if(ISI.type == "Exponential"){
      out <- x * exp(-X)
    }
    if(ISI.type == "Gamma"){
      g <- hyper
      out <- (g * x * (g*X)**(g-1) * exp(-g*X) )/(gamma(g))
    }
    if(ISI.type == "Gamma2"){
      a <- hyper[1] ; b <- hyper[2]
      out <- (b * x * (b*X)**(a-1) * exp(-b*X) )/(gamma(a))
    }
    if(ISI.type == "InverseGaussianOLD"){
      a <- hyper
      out <- (x/(sqrt(2*pi*(X**3) ))) * exp(-( (X-a)**2)/(2*X*(a**2)))
    }
    if(ISI.type == "InverseGaussian2"){
      mu <- hyper[1] ; l <- hyper[2]
      out <- (x*sqrt(l)/(sqrt(2*pi*(X**3) ))) * exp(-( l*(X-mu)**2)/(2*X*(mu**2)))
    }
    if(ISI.type == "InverseGaussian"){
      l <- hyper
      out <- x*sqrt(l/(2*pi*X**3)) * exp((-l* (X-1)**2)/(2*X))
    }
    if(ISI.type == "LogNormal"){
      mu <- hyper
      out <- (x/(2*X*sqrt(mu*pi)))*exp(-((log(X) + mu)**2)/(4*mu))
    }
    if(ISI.type == "LogNormal2"){
      mu <- hyper[1] ; sigma <- hyper[2]
      out <- (x/(sigma*X*sqrt(2*pi)))*exp(-((log(X) - mu)**2)/(2*sigma**2))
    }
    if(ISI.type == "Weibull"){
      k <- hyper
      l <- 1/(gamma(1+1/k))
      out <-   k*x/l * (X/l)**(k-1) * exp(- (X/l)**k)
    }
    if(ISI.type == "Weibull2"){
      k <- hyper[1]; l <- hyper[2]
      out <-   k*x/l * (X/l)**(k-1) * exp(- (X/l)**k)
    }
  }
  else if (do.log){
    if(ISI.type == "Exponential"){
      out <- log(x)  -X
    }
    if(ISI.type == "Gamma"){
      g <- hyper
      # out <- (g * x * (g*X)**(g-1) * exp(-g*X) )/(gamma(g))
      out <- g*log(g) + log(x) +(g-1)*log(X)  -g*X - lgamma(g)
    }
    if(ISI.type == "InverseGaussian"){
      l <- hyper
      # out <- x*sqrt(l/(2*pi*X**3)) * exp((-l* (X-1)**2)/(2*X))
      out <- log(x)+ 0.5*log(l) - 0.5*log(2*pi) - 1.5*log(X) - l* ((X-1)**2)/(2*X)
    }
    if(ISI.type == "LogNormal"){
      mu <- hyper
      out <- (x/(2*X*sqrt(mu*pi)))*exp(-((log(X) + mu)**2)/(4*mu))
      out <- log(x) - 0.5*log(2*pi) - 0.5*log(mu) - log(X) -((log(X) + mu)**2)/(4*mu)
    }
    if(ISI.type == "Weibull"){
      k <- hyper
      l <- 1/(gamma(1+1/k))
      # out <-   k*x/l * (X/l)**(k-1) * exp(- (X/l)**k)
      out <- log(x)+log(k) + (k-1)*log(X) - k*log(l) - (X/l)**k
    }
    if(ISI.type == "Gamma2"){
      a <- hyper[1] ; b <- hyper[2]
      # out <- (b * x * (b*X)**(a-1) * exp(-b*X) )/(gamma(a))
      out <- log(x) +a*log(b) +(a-1)*log(X) - b*X - lgamma(a)
    }
    if(ISI.type == "InverseGaussianOLD"){
      a <- hyper
      # out <- (x/(sqrt(2*pi*(X**3) ))) * exp(-( (X-a)**2)/(2*X*(a**2)))
      out <- log(x) - 0.5*log(2*pi) - 1.5*log(X) - (((X-mu)**2)/(2*mu**2*X))
    }
    if(ISI.type == "InverseGaussian2"){
      mu <- hyper[1] ; l <- hyper[2]
      # out <- (x*sqrt(l)/(sqrt(2*pi*(X**3) ))) * exp(-( l*(X-mu)**2)/(2*X*(mu**2)))
      out <- log(x) + 0.5*log(l) - 0.5*log(2*pi) - 1.5*log(X) - ((l*(X-mu)**2)/(2*mu**2*X))
    }
    if(ISI.type == "LogNormal2"){
      mu <- hyper[1] ; sigma <- hyper[2]
      # out <- (x/(sigma*X*sqrt(2*pi)))*exp(-((log(X) - mu)**2)/(2*sigma**2))
      out <- log(x) - log(X) - log(sigma) - 0.5*log(2*pi) - ((log(X)-u)**2)/(2*sigma**2)
    }
    if(ISI.type == "Weibull2"){
      k <- hyper[1]; l <- hyper[2]
      # out <-   k*x/l * (X/l)**(k-1) * exp(- (X/l)**k)
      out <-   log(x) + log(k) + (k-1)*log(X) - k*log(l) - (X/l)**k
    }
  }

  return(out)
}


#' Title
#'
#' @param end.time The length of time we want to simulate spikes
#' @param int.fn The intensity function
#' @param hyper The ISI parameter
#' @param steps The number of steps to split the experiment into
#' @param T.min Refractory period (Where pdf = 0)
#' @param ISI.type The ISI distribution
#' @param do.log Flag for whether we do the calculations on the log scale
#'
#' @return List where $spikes contains the spike sequence and $err gives the max error
#' @export
#'
#' @examples
Spikes <- function(end.time, int.fn, hyper, steps =1000, T.min = NULL, ISI.type = "Gamma",do.log = T){
  # Do some checks initially.
  {
    if((length(int.fn)-1)%%steps!=0){
      stop("The number of steps doesn't divide into the discretization of the intensity funtion.\n")
    }
    if(ISI.type == "Gamma" && hyper<=0){
      stop("hyper must be >0. \n")
    }
    if(max(int.fn)*(end.time/steps)>0.3){
      cat("Warning! The discretisation may not be fine enough.\n I.e max intensity =",max(int.fn),
          " and fineness = ",end.time/steps, " So expected number of spikes in interval is",
          max(int.fn)*(end.time/steps),"\n")
    }
  }

  # calculate x, X and time grid, where we step time in dt.
  dt <- end.time/steps
  x <- int.fn
  X <- cumsum(x)*(end.time/(length(x)-1))
  t <- seq(0,end.time, dt)

  # Initalise our variables.
  t <- 0
  Spikes <- 0
  last.spike <- 0
  err <- 0

  # Check for T.min (and convert to a value fo the grid if necessary).
  if(!is.null(T.min) == TRUE){
    T.min <- T.min - (T.min%%dt)
    t <- T.min
    last.spike <- T.min
  }

  while(t<(end.time-dt)){   #  We go to end.time-dt to stop having a spike at t=end.time and then entering the loop again.
    # Which would cause an eror as x is not defined for t= end.time +dt.
    t <- t + dt
    # ----------------------------
    # Loop To calcuate the next spike time (stochastic)
    xi <- stats::runif(1)
    CDF <- 0

    while(CDF<xi){
      if(CDF==0){
        f1<-0
      }else{
        if(!do.log){
          f1<- PDF(t-dt, last.spike, hyper, end.time, x, X, ISI.type,do.log)
        }else{
          f1<- exp(PDF(t-dt, last.spike, hyper, end.time, x, X, ISI.type,do.log))
        }
      }
      if(!do.log){ f2<- PDF(t, last.spike, hyper, end.time, x, X, ISI.type,do.log)}
      else{ f2<- exp(PDF(t, last.spike, hyper, end.time, x, X, ISI.type,do.log))}
      # print(f1)
      # print(f2)
      IntegralInStep <- dt*(min(f1,f2) + abs(f1-f2)/2)
      CDF<- CDF +IntegralInStep

      if(IntegralInStep>err){
        err<- IntegralInStep
      }
      if(CDF < xi){
        t<-t+dt
      }else {
        Spikes <- c(Spikes,t)
      }
      if(t>end.time){
        break
      }
    }
    last.spike <- Spikes[length(Spikes)]
    # ----------------------------

    # ----------------------------
    # If T.min is set. No spikes in the interval of length T.min after a spike.
    if (!is.null(T.min) == TRUE){
      t <- t +T.min
      last.spike <- t
    }
    # ----------------------------

  }
  return(list(error =err,spikes =Spikes[-1]))
}


#' Title
#'
#' @param end.time The length of time we want to simulate spikes
#' @param int.fn The intensity function
#' @param hyper The values of the ISI parameter
#' @param steps The number of steps to discretise time into
#' @param T.min The refractory period
#' @param ISI.type The ISI distribution
#' @param multi The number of spike sequences to generate
#' @param add.end  Do we include the end time on the end of the spikes
#' @param do.log Do we do the calculations on the log scale
#'
#' @return Matrix containing multiple spikes sequences
#' @export
#'
#' @examples
MultiSpikes <- function(end.time, int.fn, hyper, steps =1000, T.min = NULL, ISI.type = "Gamma", multi = 10, add.end = TRUE,do.log = T){
  if(add.end == TRUE){
    # Iterate over the number of sequences required.
    for(i in 1:multi){
      # Save first spike sequence in data.frame
      if(i == 1){
        V1 <- c(Spikes(end.time, int.fn, hyper, steps, T.min, ISI.type, do.log)$spikes, end.time)
        d.spikes <- data.frame(V1)
      }
      else{
        spikes <- c(Spikes(end.time, int.fn, hyper, steps, T.min, ISI.type,do.log)$spikes, end.time)

        # Add NA to end of spikes if length is less then data frame.
        if(length(spikes) < nrow(d.spikes) ){
          na.s <- rep(NA,(nrow(d.spikes)-length(spikes)))
          spikes <- c(spikes, na.s)
        }

        # Add NA to end of data frame if length spikes is larger.
        if(length(spikes) > nrow(d.spikes) ){
          d.spikes[(nrow(d.spikes)+1):(length(spikes)),] <- NA
        }
        d.spikes[,ncol(d.spikes)+1] <- spikes
      }
    }
  }
  else{
    # Iterate over the number of sequences required.
    for(i in 1:multi){
      # Save first spike sequence in data.frame
      if(i == 1){
        V1 <- Spikes(end.time, int.fn, hyper, steps, T.min, ISI.type,do.log)$spikes
        d.spikes <- data.frame(V1)
      }
      else{
        spikes <- Spikes(end.time, int.fn, hyper, steps, T.min, ISI.type,do.log)$spikes

        # Add NA to end of spikes if length is less then data frame.
        if(length(spikes) < nrow(d.spikes) ){
          na.s <- rep(NA,(nrow(d.spikes)-length(spikes)))
          spikes <- c(spikes, na.s)
        }

        # Add NA to end of data frame if length spikes is larger.
        if(length(spikes) > nrow(d.spikes) ){
          d.spikes[(nrow(d.spikes)+1):(length(spikes)),] <- NA
        }
        d.spikes[,ncol(d.spikes)+1] <- spikes
      }
    }
  }

  return(d.spikes)
}
