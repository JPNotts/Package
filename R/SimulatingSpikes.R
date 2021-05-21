
#' Simulate a single spike sequence
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
#'
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


#' Simulate spike sequences
#'
#' @inheritParams Spikes
#' @param sequences The number of spike sequences to generate
#' @param add.end  Do we include the end time on the end of the spikes
#'
#' @return Matrix containing multiple spikes sequences
#' @export
#'
#' @examples
simulate_spikes <- function(end.time, int.fn, hyper, steps =1000, T.min = NULL, ISI.type = "Gamma", sequences = 1, add.end = TRUE,do.log = T){
  if(add.end == TRUE){
    # Iterate over the number of sequences required.
    for(i in 1:sequences){
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
    for(i in 1:sequences){
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
