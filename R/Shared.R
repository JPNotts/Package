#' Calculate the probability density.
#'
#' @param t The time of the current spike
#' @param last.spike The time of the previous spike.
#' @param hyper The ISI parameter value
#' @param end.time The end time of the experiment
#' @param x The intensity function --- defined in 0,end.time
#' @param X The integral of the intensity function
#' @param ISI.type The ISI distribution.
#' @param do.log Flag for whether to calculate the pdf on the log scale.
#'
#' @return The probability density of a spike at time t given the last spike was at last.spike.
#' @export
#'
#' @examples
PDF <- function(t, last.spike, hyper, end.time, x, X, ISI.type, do.log = F){
  # Check the ISI type is allowable
  ISIs <- c('Gamma', 'Exponential', 'InverseGaussian', 'LogNormal', 'Weibull',
            'InverseGaussianOLD','Gamma2', 'InverseGaussian2', 'LogNormal2', 'Weibull2' )
  if(!(ISI.type %in% ISIs)){ stop("Invalid ISI.type. See help documentation for allowable ISI.type.")}
  ISI2 <- c('Gamma2', 'InverseGaussian2', 'LogNormal2', 'Weibull2')
  if(ISI.type %in% ISI2 && length(hyper) != 2){stop("Your choice of ISI.type requires two inputs for hyper.")}
  ISI2 <- c('Gamma', 'InverseGaussian', 'LogNormal', 'Weibull','InverseGaussianOLD')
  if(ISI.type %in% ISI2 && length(hyper) != 1){stop("Your choice of ISI.type requires one inputs for hyper. For example hyper=1.")}

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
      out <- log(x) - log(X) - log(sigma) - 0.5*log(2*pi) - ((log(X)-mu)**2)/(2*sigma**2)
    }
    if(ISI.type == "Weibull2"){
      k <- hyper[1]; l <- hyper[2]
      # out <-   k*x/l * (X/l)**(k-1) * exp(- (X/l)**k)
      out <-   log(x) + log(k) + (k-1)*log(X) - k*log(l) - (X/l)**k
    }
  }

  return(out)
}


#' Calculate the probability density.
#'
#'
#'
#' @param t The time of the current spike
#' @param last.spike The time of the previous spike.
#' @param hyper The ISI parameter value
#' @param end.time The end time of the experiment
#' @param x The intensity function --- defined in 0,end.time
#' @param X The integral of the intensity function
#' @param ISI.type The ISI distribution.
#' @param do.log Flag for whether to calculate the pdf on the log scale.
#' @param t.min t.min
#'
#' @return The probability density of a spike at time t given the last spike was at last.spike.
#' @export
#'
#' @examples
PDF_tmin <- function(t, last.spike, hyper, end.time, x, X, ISI.type, t.min = 0, do.log = F){
  # Check the ISI type is allowable
  ISIs <- c('Gamma', 'Exponential', 'InverseGaussian', 'LogNormal', 'Weibull',
            'InverseGaussianOLD','Gamma2', 'InverseGaussian2', 'LogNormal2', 'Weibull2' )
  if(!(ISI.type %in% ISIs)){ stop("Invalid ISI.type. See help documentation for allowable ISI.type.")}
  ISI2 <- c('Gamma2', 'InverseGaussian2', 'LogNormal2', 'Weibull2')
  if(ISI.type %in% ISI2 && length(hyper) != 2){stop("Your choice of ISI.type requires two inputs for hyper.")}
  ISI2 <- c('Gamma', 'InverseGaussian', 'LogNormal', 'Weibull','InverseGaussianOLD')
  if(ISI.type %in% ISI2 && length(hyper) != 1){stop("Your choice of ISI.type requires one inputs for hyper. For example hyper=1.")}



  # First if t-last.spike < t.min return 0
  if(t-last.spike < t.min){
    return(0)
  }
  step.size <- end.time/(length(x)-1)
  x <- x[t/step.size]
  if(last.spike - 0 < 1e-10){
    
    X <- X[t/step.size]
  }
  else{
    X <- X[t/step.size] - X[(last.spike + t.min)/step.size]
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
      out <- log(x) - log(X) - log(sigma) - 0.5*log(2*pi) - ((log(X)-mu)**2)/(2*sigma**2)
    }
    if(ISI.type == "Weibull2"){
      k <- hyper[1]; l <- hyper[2]
      # out <-   k*x/l * (X/l)**(k-1) * exp(- (X/l)**k)
      out <-   log(x) + log(k) + (k-1)*log(X) - k*log(l) - (X/l)**k
    }
  }

  return(out)
}

#' get_spikes()
#'
#'A function that returns a single spike sequence from a data frame of spikes.
#' @param d.spikes Data frame containing spike sequences.
#' @param i  The index (column) denoting the required spike sequence.
#' @param rm.end.time  A flag, where if true we remove the last spike from the spike sequence --- encase it corresponds to the experiment end time.
#'
#' @return A single spike sequence in a vector.
#' @export
#'
#'
get_spikes <- function(d.spikes,i, rm.end.time = TRUE){
  # First check i is feasible.
  if(i > ncol(d.spikes)){
    stop("You are trying to get the", i,"th spike sequence but you only have", ncol(d.spikes),"spike sequences!")
  }
  spikes <- d.spikes[,i][!is.na(d.spikes[,i])]
  if(rm.end.time == TRUE){
    spikes <- spikes[-length(spikes)]
  }

  return(spikes)
}

#' get_X()
#'
#'Calculates the integral of the intensity function, where the intensity is a piece-wise constant function.
#'
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param d.spikes Data frame containing spike sequences
#' @param T.min Refractory period, default set to NULL.
#' @param rm.end.time Flag, where if TRUE we remove the end time of the spike sequences.
#' @param x The intensity function
#' @param end.time End of experiment time
#' @param int.method  Method for integration
#'
#' @return data frame containing the integral of the intensity between all spike times.
#' @export
#'
get_X <- function(d.spikes,partition=NA, heights=NA, x=NA, end.time=NA, T.min = NULL, rm.end.time=T, int.method = "trapezium")
{
  # GP/Constant get_X.
  if(!is.na(x[1])){

    # If x is constant (i.e length =1)
    if(length(x) == 1){
      X <- matrix(0, ncol=dim(d.spikes)[2], nrow = (dim(d.spikes)[1]))
      intensity <- matrix(NA, ncol=dim(d.spikes)[2], nrow = (dim(d.spikes)[1]-1))

      for(m in 1:ncol(d.spikes)){
        # Get the current spike sequence.
        spikes <- get_spikes(d.spikes,m)

        # Due to allowing multiple spike sequences they may vary in length. Calculate diff to concatenate NA to end of X and intensity.
        s <- length(spikes)
        diff <- length(intensity[,1]) - s

        # First compute the intensity
        intensity[,m] <- c(rep(x,s), rep(NA,diff))

        # Next we want to get X.
        pre.spikes <- spikes
        if(!is.null(T.min)){pre.spikes <- spikes +T.min}
        pre <-c(0,pre.spikes) ; post <-c(spikes,end.time)
        X[,m] <-c( x*(post - pre), rep(NA,diff))

        # Check if the last X entry is negative (when end of experiment occurs before the end of refractory period).
        # If it is set to NA.
        if(X[length(spikes),m] <0 ){
          X[length(spikes),m] <- NA
        }
      }
      return(list(integral = X, intensity = intensity))

    }

      # First we need to create a store for X and x allowing for multiple spike sequences.
      # Choose to store in matrices where each column corresponds to the values for the spike sequence d.spikes[,i].
      X <- matrix(0, ncol=dim(d.spikes)[2], nrow = (dim(d.spikes)[1]))
      intensity <- matrix(NA, ncol=dim(d.spikes)[2], nrow = (dim(d.spikes)[1]-1))

      # Create X.vect which is the integral of x.
      if(int.method == "trapezium"){
        add <- (x[2:length(x)] +x[1:(length(x)-1)])/2
        X.vect <- c(0,cumsum(add)*(end.time/(length(x)-1)))
      }
      # Else if we choose to use box style
      else if(int.method == "box"){
        X.vect <- c(0,cumsum(x)*(end.time/(length(x)-1)))
      }
      # else print error message
      else{
        stop("You haven't entered a correct type of integral.")
      }

      # Loop over all the spike sequences.
      for(m in 1:ncol(d.spikes)){
        # Get the current spike sequence.
        spikes <- get_spikes(d.spikes,m)

        # Due to allowing multiple spike sequences they may vary in length. Caluclate diff to concaternate NA to end of X and intensity.
        diff <- length(intensity[,1]) - length(spikes)

        # First compute the intensity.
        position <- spikes/(end.time/(length(x)-1)) +1
        intensity[,m] <- c(x[position], rep(NA,diff))

        # Next we want to get X.
        pre.spikes <- spikes
        if(!is.null(T.min)){pre.spikes <- spikes +T.min}
        pre <-c(0,pre.spikes)/(end.time/(length(x)-1))+1 ; post <-c(spikes,end.time)/(end.time/(length(x)-1))+1
        X[,m] <-c(X.vect[post] - X.vect[pre], rep(NA,diff))

        # Check if the last X entry is negative (when end of experiment occurs before the end of refractoy period).
        # If it is set to NA.
        if(X[length(spikes),m] <0 ){
          X[length(spikes),m] <- NA
        }
      }
      return(list(integral = X, intensity = intensity))
  }

  # PWC get_X.
  else if(!is.na(partition[1]) && !is.na(heights[1])){
    # First we need to create a store for X and x allowing for multiple spike sequences.
    # Choose to store in matrices where each column corresponds to the values for the spike sequence d.spikes[,i].

    # We first need to be careful, and check whether the end time is also in the spikes.
    # If it is we need to nake nrow = dim(d.spikes)[1] rather than nrow = dim(d.spikes)[1] +1
    if(abs(partition[length(partition)] - max(d.spikes[!is.na(d.spikes)])) < 1e-10){
      rows <- dim(d.spikes)[1]
      rm.end.time <- TRUE
    } else{rows = dim(d.spikes)[1] +1
    rm.end.time = F
    }

    X <- matrix(0, ncol=dim(d.spikes)[2], nrow = rows)
    intensity <- matrix(NA, ncol=dim(d.spikes)[2], nrow = (rows-1))

    # Loop over all the spike sequences.
    for(m in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      spikes <- get_spikes(d.spikes,m,rm.end.time)

      # Compute X and x if T.min is NULL.
      if(is.null(T.min) == TRUE){
        together <- sort(c(partition,spikes))

        # Counters, j changes the integral we're dealing with, k changes the height value.
        j<-1
        k <- 1

        h.cur <- heights[k]
        for(i in 2:length(together)){
          if(together[i] == together[i-1]){ next }
          X.in.step <- h.cur * (together[i]-together[i-1])
          X[j,m] <- X[j,m] + X.in.step

          if (together[i] %in% partition){
            k <- k+1
            h.cur <- heights[k]
          }
          if (together[i] %in% spikes){
            intensity[j,m] <- h.cur
            j <- j+1
          }
        }

      }
      else{ # Compute X and x if we have a T.min value.

        # Calculate the values of last spike + T.min.
        spikes.with.T.min <- spikes + T.min

        # If the experiment ends in the T.min period after the last spike set X(N+1) = NA.
        if (spikes.with.T.min[length(spikes)] > partition[length(partition)]){
          X[length(X),m] <- NA
          spikes.with.T.min <- spikes.with.T.min[-length(spikes)]
        }

        together <- sort(c(partition,spikes,spikes.with.T.min))

        # Counters, j changes the integral we're dealing with, k changes the height value.
        j<-1
        k <- 1
        allow <- TRUE  # A FLAG that tells you whether to calculate and add the integral in the step.
        h.cur <- heights[k]

        for(i in 2:length(together)){
          if(together[i] == together[i-1]){ next }
          if(allow == TRUE){
            X.in.step <- h.cur * (together[i]-together[i-1])
            X[j,m] <- X[j,m] + X.in.step
          }

          if (together[i] %in% partition){
            k <- k+1
            h.cur <- heights[k]
          }
          if (together[i] %in% spikes){
            intensity[j,m] <- h.cur
            allow <- FALSE
            j <- j+1
          }
          if(together[i] %in% spikes.with.T.min){
            allow <- TRUE
          }
        }
      }
    }
    return(list(integral = X, intensity = intensity))
  }
  else{
    stop("You haven't entered a valid intensity function.")}
}

# A function to check whether the partition, heights and spikes are all consistent.
#' Check a step function is feasible.
#'
#'A function that checks that the partition and heights of a step function is feasible to be used as an intensity function. The function produces an error if the intensity is not feasible.
#'
#' @inheritParams get_X
#' @param max.T.min  If T.min != NULL, the maximum allowable refractory period.
#'
#' @export
check_feasible <- function(d.spikes, partition=NA ,heights=NA, x=NA, end.time = NA, T.min = NULL, max.T.min=NA){
  # Check partition/heights based issues.
  if(!is.na(partition[1]) && !is.na(heights[1]) && is.na(x[1]) && is.na(end.time)){

      # Firstly check that partition and heights has the correct form.
      if(length(partition) != length(heights) + 1){
        stop("Error! Your partition of time has ",length(partition)," values and you have", length(heights),
             " height values, whereas you should have ",length(partition)-1," values. /n")
      }

      # Check that partition is in increasing order.
      if ( any(partition != sort(partition))){
        stop("Error! The partition is not in ascending order.\n partition = ",partition,"\n")
      }

  }
  # Check x/end.time based issues.
  else if (is.na(partition[1]) && is.na(heights[1]) && !is.na(x[1]) && !is.na(end.time)){
    # Check x is always positive
    if(any(x < 0)){
      stop("Intensity contains negative values!")
    }
    # Check end.time is after the last spike time.
    if(max(unlist(d.spikes),na.rm=T) > end.time +1e-10){
      stop("Spikes outside experiment time!")
    }
  }
  else{
    stop('Need values for only partition and heights or x and end.time!')
  }
  # Check t.min value is allowable.
  if(!is.null(T.min)){
    # Make sure max.T.min has a value
    if(  !(length(max.T.min) > 0 && is.numeric(max.T.min)) ){
      stop("max.T.min value not valid!")
    }
    # Check T.min not greater than max.T.min
    if(T.min > max.T.min){
      stop("T.min is larger then maximal allowed value!")
    }
  }

}


#' Likelihood of the intensity function.
#'
#' @inheritParams check_feasible
#' @param hyper.param Value of the ISI parameter.
#' @param ISI.type The ISI distribution.
#' @param do.log Flag, where if do.log is true the calculations are computed on the log scale.
#'
#' @return Likelihood of the intensity function.
#' @export
#'
log_pi_x <- function(d.spikes, hyper.param, partition = NA, heights= NA, x= NA, end.time = NA, T.min = NULL, max.T.min = NA, ISI.type = "Gamma", do.log = TRUE){

  # Intensity from PWC
  if(!is.na(partition[1]) && !is.na(heights[1]) && is.na(x[1]) && is.na(end.time)){
    # Firstly check everything is consistent.
    check_feasible(d.spikes,partition = partition, heights = heights, T.min = T.min, max.T.min = max.T.min)
    values <- get_X(d.spikes,partition = partition, heights = heights, T.min = T.min)
    d.x <- values$intensity  ; d.X <- values$integral
  }
  # Intensity from GP
  else if (is.na(partition[1]) && is.na(heights[1]) && !is.na(x[1]) && !is.na(end.time)){
    check_feasible(d.spikes, x = x, end.time = end.time, T.min = T.min, max.T.min = max.T.min)
    values <- get_X(d.spikes, x=x, end.time=end.time, T.min=T.min)
    d.x <- values$intensity  ; d.X <- values$integral
  }
  else{stop('Require valid choice between GP and PWC.')}

  ##----------------------------
  ## FOR IN GAMMA ISI DISTRIBUTION.
  ##----------------------------
  if(ISI.type == "Gamma"){
    gamma <- hyper.param

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.

      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)
      # print('N')
      # print(N)
      # print('x')
      # print(x)
      # print('X')
      # print(X)
      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          # print('boonah :(')
          out <- out + log(x[1]) - X[1]
        }
        else{
          # print('booyah!')
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(gamma) + base::log(x[2:N]) - lgamma(gamma)  + (gamma-1)*base::log(gamma*X[2:N]) - gamma * X[2:N])

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( ( gamma * x[2:N] / gamma(gamma) ) * ((gamma*X[2:N] )**(gamma-1)) * exp(-gamma * X[2:N]) )
      }
    }

  }
  ##----------------------------
  ## FOR INHOMO INVERSE GAUSSIAN DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "InvGauss"){
    alpha <- hyper.param

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(x[2:N]) -0.5*base::log(2*pi) -1.5*base::log(X[2:N]) - sum((X[2:N] - alpha)**2)/(2*X[2:N]*alpha**2))
      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( (x[2:N]/ ((2*pi * (X[2:N]**3) )**0.5) ) * exp(-((X[2:N] - alpha)**2)/(2*(alpha**2) * X[2:N])) )
      }
    }
  }
  ##----------------------------
  ## FOR New INVERSE GAUSSIAN DISTRIBUTION (mean 1).
  ##----------------------------
  else if(ISI.type == "NewIG"){
    l <- hyper.param

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(x[2:N]) + 0.5*base::log(l) - 0.5*base::log(2*pi) - 1.5*base::log(X[2:N]) - (l*(X[2:N]-1)**2)/(2*X[2:N]))

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( ((x[2:N]*sqrt(l))/ ((2*pi * (X[2:N]**3) )**0.5) ) * exp(-l*((X[2:N] - 1)**2)/(2* X[2:N])) )
      }
    }
  }
  ##----------------------------
  ## FOR INHOMO POISSON.
  ##----------------------------
  else if(ISI.type == "Poisson"){

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.
    if(do.log == TRUE){out <- 0}else{out <- 1}

    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(x[2:N]) - X[2:N])

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( x[2:N] * exp(-X[2:N]) )
      }
    }
  }
  ##----------------------------
  ##----------------------------
  ## FOR LOG NORMAL DISTRIBUTION (mean 1).
  ##----------------------------
  else if(ISI.type == "LN"){
    mu <- hyper.param

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(x[2:N]) -base::log(2*X[2:N]) - 0.5*base::log(mu*pi) - (base::log(X[2:N]) +mu)**2/(4*mu))

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( (x[2:N]/ (2*X[2:N] *sqrt(mu*pi)) ) * exp( -((log(X[2:N])+ mu)**2)/(4*mu) ))
      }
    }

  }
  ##----------------------------
  ##----------------------------
  ## FOR WEIBULL DISTRIBUTION (mean 1).
  ##----------------------------
  else if(ISI.type == "Weibull"){
    k <- hyper.param
    l <- 1/(gamma(1+1/k))
    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous Weibull distribution)
        out <- out + sum( base::log(x[2:N]) + base::log(k/l) +(k-1)*base::log(X[2:N]/l) - (X[2:N]/l)**k)

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod((k*x[2:N]/l)*(X[2:N]/l)**(k-1) * exp(-(X[2:N]/l)**k) )
      }
    }
  }
  ##----------------------------
  else{
    stop("Error! You haven't chose a correct ISI.type. \n You inputted:\"", ISI.type,"\".\n")}


  return(out)
}


#' Marginal of the ISI parameter
#'
#' @inheritParams log_pi_x
#' @param prior.hyper.param Prior of the hyper parameter.
#'
#' @return marginal of the ISI parameter
#' @export
log_pi_hyper_param <- function(d.spikes, hyper.param, partition=NA, heights=NA, x=NA, end.time = NA, prior.hyper.param = c(1,0.01), T.min = NULL, max.T.min = NA, ISI.type = "Gamma"){
  # Intensity from PWC
  if(!is.na(partition[1]) && !is.na(heights[1]) && is.na(x[1]) && is.na(end.time)){
    # Firstly check everything is consistent.
    check_feasible(d.spikes,partition = partition, heights = heights, T.min = T.min, max.T.min = max.T.min)
    values <- get_X(d.spikes,partition = partition, heights = heights, T.min = T.min)
    d.x <- values$intensity  ; d.X <- values$integral
  }
  # Intensity from GP
  else if (is.na(partition[1]) && is.na(heights[1]) && !is.na(x[1]) && !is.na(end.time)){
    check_feasible(d.spikes, x = x, end.time = end.time, T.min = T.min, max.T.min = max.T.min)
    values <- get_X(d.spikes, x=x, end.time=end.time, T.min=T.min)
    d.x <- values$intensity  ; d.X <- values$integral
  }


  # Initially add the value from the prior
  out <- (prior.hyper.param[1]-1) * base::log(hyper.param) - prior.hyper.param[2]*hyper.param

  ##----------------------------
  ## FOR IN GAMMA ISI DISTRIBUTION.
  ##----------------------------
  if(ISI.type == "Gamma"){
    gam <- hyper.param
    if (gam < 0){
      stop("Error! Your gamma value is less then 0. gamma = ",gam,".")
    }

    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i, rm.end.time = F)
      N <- length(X) - 1

      out <- out +(N-1)*( gam * base::log(gam) - lgamma(gam))
      out <- out + sum( (gam-1)* base::log(X[2:N]) - gam*X[2:N] )
    }
  }
  ##----------------------------
  ## FOR IN INVERSE GAUSSIAN ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "InvGauss"){
    alpha <- hyper.param
    if (alpha < 0){
      stop("Error! Your alpha value is less then 0. alpha = ",alpha,".")
    }

    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i, rm.end.time = F)
      x <- get_spikes(d.x,i, rm.end.time = FALSE)
      N <- length(X) - 1
      out <- out +  sum( base::log(x[2:N]) - 0.5*base::log(2*pi) -3/2*base::log(X[2:N])  -( (X[2:N]-alpha)**2)/(2*X[2:N]*(alpha**2)) )
      # out <- out + sum( (-(X[2:N]-alpha)**2)/(2*alpha**2 * X[2:N]) )
    }


  }
  ##----------------------------
  ## FOR IN INVERSE GAUSSIAN ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "NewIG"){
    l <- hyper.param
    if (l < 0){
      stop("Error! Your l value is less then 0. l = ",l,".")
    }

    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i, rm.end.time = F)
      N <- length(X) - 1

      out <- out + (N-1)*base::log(l)/2 + sum( (-l*(X[2:N]-1)**2)/(2*X[2:N]) )
    }
  }
  ##----------------------------
  ## FOR LOG NORMAL ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "LN"){
    mu <- hyper.param
    if (mu < 0){
      stop("Error! Your mu value is less then 0. mu = ",mu,".")
    }
    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i, rm.end.time = FALSE)
      x <- get_spikes(d.x,i, rm.end.time = FALSE)
      N <- length(X) - 1
      out <- out + sum(base::log(x[2:N]/(2*X[2:N]*sqrt(pi*mu))) - ((base::log(X[2:N]) +mu)**2)/(4*mu) )
    }
  }
  ##----------------------------
  ## FOR WEIBULL ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "Weibull"){
    k <- hyper.param
    if (k < 0){
      stop("Error! Your k value is less then 0. k = ",k,".")
    }
    if (k< 0.06){
      out <- -Inf
    }
    else{
      l <- 1/(gamma(1+1/k))

      # Iterate over all spike sequences adding each contribution.
      for(i in 1:ncol(d.spikes)){
        X <- get_spikes(d.X,i, rm.end.time = F)
        N <- length(X) - 1

        out <- out + sum( base::log(k/l) +(k-1)*base::log(X[2:N]/l) - (X[2:N]/l)**k)
        # out <- out + (N-1)*(k*lgamma(1+1/k) +log(k)) +sum( (k-1)*log(X[2:N]) - (gamma(1+1/k)*X[2:N])**k)
      }


    }
  }
  ##----------------------------
  else{
    stop("Error! You haven't chose a correct ISI.type. \n You inputted:\"", ISI.type,"\".\n")
  }

  return(out)
}


#' Likelihood of the refractory period
#'
#' @inheritParams log_pi_x
#' @param prior.T.min Prior for the refractory period.
#'
#' @return marginal of the refractory period
#' @export
#'
log_pi_Tmin <- function(d.spikes, T.min, max.T.min, hyper.param, partition = NA, heights = NA, x=NA, end.time=NA, prior.T.min = c(1,0.01), ISI.type = "Gamma"){
  # Intensity from PWC
  if(!is.na(partition[1]) && !is.na(heights[1]) && is.na(x[1]) && is.na(end.time)){
    # Firstly check everything is consistent.
    check_feasible(d.spikes,partition = partition, heights = heights, T.min = T.min, max.T.min = max.T.min)
    values <- get_X(d.spikes,partition = partition, heights = heights, T.min = T.min)
    d.x <- values$intensity  ; d.X <- values$integral
  }
  # Intensity from GP
  else if (is.na(partition[1]) && is.na(heights[1]) && !is.na(x[1]) && !is.na(end.time)){
    check_feasible(d.spikes, x= x, end.time = end.time, T.min = T.min, max.T.min = max.T.min)
    values <- get_X(d.spikes, x=x, end.time=end.time, T.min=T.min)
    d.x <- values$intensity  ; d.X <- values$integral
  }


  # calculate the prior contribution.
  out <- (prior.T.min[1]-1) * log(T.min) - prior.T.min[2]*T.min

  # Next add the contribution for p1 and pT
  for(i in 1:ncol(d.spikes)){
    X <- get_spikes(d.X,i)
    if(is.na(X[length(X)])){ out <- out - X[1]}
    else{ out <- out - X[1] - X[length(X)]}
  }


  ##----------------------------
  ## FOR POISSON ISI DISTRIBUTION.
  ##----------------------------
  if(ISI.type == "Poisson"){
    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i)  # although the function is get spikes this returns the X[] for the sequence.
      N <- length(X) - 1

      out <- out + sum( - X[2:N] )
    }
  }
  ##----------------------------
  ##----------------------------
  ## FOR IN GAMMA ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "Gamma"){
    gam <- hyper.param
    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i)  # although the function is get spikes this returns the X[] for the sequence.
      N <- length(X) - 1

      out <- out + sum( (gam-1)* log(X[2:N]) - gam*X[2:N] )
    }
  }
  ##----------------------------
  ## FOR IN INVERSE GAUSSIAN ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "InverseGaussian"){
    alpha <- hyper.param
    if (alpha < 0){
      stop("Error! Your alpha value is less then 0. alpha = ",alpha,".")
    }
    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i)
      N <- length(X) - 1

      out <- out + sum( (-(X[2:N]-alpha)**2)/(2*alpha**2 * X[2:N]) )
    }
  }
  ##----------------------------
  ## FOR IN INVERSE GAUSSIAN ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "NewIG"){
    l <- hyper.param
    if (l < 0){
      stop("Error! Your l value is less then 0. l = ",l,".")
    }
    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i)
      N <- length(X) - 1

      out <- out + sum( (-l*(X[2:N]-1)**2)/(2*X[2:N]) )
    }
  }
  ##----------------------------
  ## FOR LOG NORMAL ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "LN"){
    mu <- hyper.param
    if (mu < 0){
      stop("Error! Your mu value is less then 0. mu = ",mu,".")
    }
    d.x <- get_X(partition,heights,d.spikes,T.min)$intensity

    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i, rm.end.time = FALSE)
      x <- get_spikes(d.x,i, rm.end.time = FALSE)
      N <- length(X) - 1
      out <- out + sum(log(x[2:N]/(2*X[2:N]*sqrt(pi*mu))) - ((log(X[2:N]) +mu)**2)/(4*mu) )
    }
  }
  ##----------------------------
  ## FOR WEIBULL ISI DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "Weibull"){
    k <- hyper.param
    if (k < 0){
      stop("Error! Your k value is less then 0. k = ",k,".")
    }
    l <- 1/(gamma(1+1/k))

    # Iterate over all spike sequences adding each contribution.
    for(i in 1:ncol(d.spikes)){
      X <- get_spikes(d.X,i)
      N <- length(X) - 1

      out <- out + sum((k-1)*log(X[2:N]/l) - (X[2:N]/l)**k)
    }
  }
  ##----------------------------
  else{
    stop("Error! You haven't chose a correct ISI.type. \n You inputted:\"", ISI.type,"\".\n")
  }

  return(out)
}


#' Get full likelihood.
#'@inheritParams log_pi_x
#'
#' @return
#' @export
#'
#' @examples
log_likelihood <- function(d.spikes, hyper.param, partition = NA, heights= NA, x= NA, end.time = NA, T.min = NULL, max.T.min = NA, ISI.type = "Gamma", do.log = TRUE){

  # Intensity from PWC
  if(!is.na(partition[1]) && !is.na(heights[1]) && is.na(x[1]) && is.na(end.time)){
    # Firstly check everything is consistent.
    check_feasible(d.spikes,partition = partition, heights = heights, T.min = T.min, max.T.min = max.T.min)
    values <- get_X(d.spikes,partition = partition, heights = heights, T.min = T.min)
    d.x <- values$intensity  ; d.X <- values$integral
  }
  # Intensity from GP
  else if (is.na(partition[1]) && is.na(heights[1]) && !is.na(x[1]) && !is.na(end.time)){
    check_feasible(d.spikes, x = x, end.time = end.time, T.min = T.min, max.T.min = max.T.min)
    values <- get_X(d.spikes, x=x, end.time=end.time, T.min=T.min)
    d.x <- values$intensity  ; d.X <- values$integral
  }
  else{stop('Require valid choice between GP and PWC.')}

  ##----------------------------
  ## FOR IN GAMMA ISI DISTRIBUTION.
  ##----------------------------
  if(ISI.type == "Gamma"){
    gamma <- hyper.param

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.

      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)
      # print('N')
      # print(N)
      # print('x')
      # print(x)
      # print('X')
      # print(X)
      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          # print('boonah :(')
          out <- out + log(x[1]) - X[1]
        }
        else{
          # print('booyah!')
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(gamma) + base::log(x[2:N]) - lgamma(gamma)  + (gamma-1)*base::log(gamma*X[2:N]) - gamma * X[2:N])

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( ( gamma * x[2:N] / gamma(gamma) ) * ((gamma*X[2:N] )**(gamma-1)) * exp(-gamma * X[2:N]) )
      }
    }

  }
  ##----------------------------
  ## FOR INHOMO INVERSE GAUSSIAN DISTRIBUTION.
  ##----------------------------
  else if(ISI.type == "InvGauss"){
    alpha <- hyper.param

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(x[2:N]) -0.5*base::log(2*pi) -1.5*base::log(X[2:N]) - sum((X[2:N] - alpha)**2)/(2*X[2:N]*alpha**2))
      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( (x[2:N]/ ((2*pi * (X[2:N]**3) )**0.5) ) * exp(-((X[2:N] - alpha)**2)/(2*(alpha**2) * X[2:N])) )
      }
    }
  }
  ##----------------------------
  ## FOR New INVERSE GAUSSIAN DISTRIBUTION (mean 1).
  ##----------------------------
  else if(ISI.type == "NewIG"){
    l <- hyper.param

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(x[2:N]) + 0.5*base::log(l) - 0.5*base::log(2*pi) - 1.5*base::log(X[2:N]) - (l*(X[2:N]-1)**2)/(2*X[2:N]))

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( ((x[2:N]*sqrt(l))/ ((2*pi * (X[2:N]**3) )**0.5) ) * exp(-l*((X[2:N] - 1)**2)/(2* X[2:N])) )
      }
    }
  }
  ##----------------------------
  ## FOR INHOMO POISSON.
  ##----------------------------
  else if(ISI.type == "Poisson"){

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.
    if(do.log == TRUE){out <- 0}else{out <- 1}

    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(x[2:N]) - X[2:N])

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( x[2:N] * exp(-X[2:N]) )
      }
    }
  }
  ##----------------------------
  ##----------------------------
  ## FOR LOG NORMAL DISTRIBUTION (mean 1).
  ##----------------------------
  else if(ISI.type == "LN"){
    mu <- hyper.param

    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out + sum(base::log(x[2:N]) -base::log(2*X[2:N]) - 0.5*base::log(mu*pi) - (base::log(X[2:N]) +mu)**2/(4*mu))

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod( (x[2:N]/ (2*X[2:N] *sqrt(mu*pi)) ) * exp( -((log(X[2:N])+ mu)**2)/(4*mu) ))
      }
    }

  }
  ##----------------------------
  ##----------------------------
  ## FOR WEIBULL DISTRIBUTION (mean 1).
  ##----------------------------
  else if(ISI.type == "Weibull"){
    k <- hyper.param
    l <- 1/(gamma(1+1/k))
    # Set out which stores our likelihood value. If we do log we add thigns together so set to 0,
    # if we don't log set to 1 as we multiply things together.

    if(do.log == TRUE){out <- 0}else{out <- 1}
    # Iterate over the spike sequences you have inputted.
    for(i in 1:ncol(d.spikes)){
      # Get the current spike sequence.
      x <- get_spikes(d.x, i, rm.end.time = F)
      X <- get_spikes(d.X, i, rm.end.time = F)

      N <- length(x)

      # Decide whether to calculate the log or not?
      if(do.log == TRUE){
        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out + base::log(x[1]) - X[1]
        }
        else{
          out <- out + base::log(x[1]) - X[1] - X[N+1]
        }

        #  calcuate all the intermediate points (inhonogeneous Weibull distribution)
        out <- out + sum( base::log(x[2:N]) + base::log(k/l) +(k-1)*base::log(X[2:N]/l) - (X[2:N]/l)**k)

      }
      else{  # Don't compute the log.

        # calculate probability for the first spike and interval after last spike.
        # Where we check to see if the interval after the last spike is important(i.e T-y_n > T.min), and
        # calculate probability for the first spike and interval after last spike.
        if(is.na(X[N+1]) == TRUE){
          out <- out * x[1] * exp(-X[1])
        }
        else{
          out <- out * x[1] * exp(-X[1]) * exp(-X[N+1])
        }

        #  calcuate all the intermediate points (inhonogeneous gamma distribution)
        out <- out * prod((k*x[2:N]/l)*(X[2:N]/l)**(k-1) * exp(-(X[2:N]/l)**k) )
      }
    }
  }
  ##----------------------------
  else{
    stop("Error! You haven't chose a correct ISI.type. \n You inputted:\"", ISI.type,"\".\n")}


  return(out)
}

