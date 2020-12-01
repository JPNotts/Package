
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

# A function to check whether the partition, heights and spikes are all consistent.
#' Check a step function is feasible.
#'
#'A function that checks that the partition and heights of a step function is feasible to be used as an intensity function. The function produces an error if the intensity is not feasible.
#'
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param d.spikes Data frame containing the spikes.
#' @param T.min Refractory period, default set to NULL.
#'
#'
#' @export
check_feasible <- function(partition,heights,d.spikes,T.min = NULL){
  # Firstly check things not to do with the spikes.
  {
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
  # Loop through each spike sequence and check the last spike is within the partition and ISI>T.min for all ISI.
  for(i in 1:ncol(d.spikes)){
    spikes <- get_spikes(d.spikes, i)

    # Check that the last spike time is within the partition.
    if(spikes[length(spikes)] > partition[length(partition)]){
      stop("Error! The last spike time is outside of the intensity function's definition.\n i.e last spike =",
           spikes[length(spikes)], " and intensity function defined upto ",partition[length(partition)], ".\n")
    }

    # If T.min != NULL, check that each ISI is at least T.min in length.
    if(!is.null(T.min) == TRUE){
      N <- length(spikes)
      min.gap <- min(spikes[2:N] - spikes[1:(N-1)])
      if( min.gap < T.min){
        stop("Error! T.min is too large! \n I.e the smallest ISI is ", min.gap," and you have set T.min to be ", T.min,". \n ")
      }
    }
  }

}


#' get_X()
#'
#'Calculates the integral of the intensity function, where the intensity is a piece-wise constant function.
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param d.spikes Data frame containing spike sequences
#' @param T.min Refractory period, default set to NULL.
#' @param rm.end.time Flag, where if TRUE we remove the end time of the spike sequences.
#'
#' @return data frame containing the integral of the intensity between all spike times.
#' @export
#'
get_X <- function(partition, heights, d.spikes, T.min = NULL, rm.end.time=T){
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


#' Poisson conditioned on max value.
#'
#' @param k k
#' @param lambda Parameter of the Poisson distribution.
#' @param k.max Maximum number of change points in RJMCMC.
#' @param sum.of.prob sum^k.max_i=0 Poi(k | lambda).
#'
#' @return probability of k
#' @export
#'
conditioned_poisson<- function(k, lambda, k.max, sum.of.prob){
  if(k > k.max){
    prob <- 0
  }
  else{
    prob <- stats::dpois(k,lambda) / sum.of.prob
  }
  return(prob)
}


#' Likelihood of the intensity function.
#'
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param d.spikes Data frame containing spike sequences.
#' @param hyper.param Value of the ISI parameter.
#' @param T.min Refractory period, default set to NULL.
#' @param ISI.type The ISI distribution.
#' @param do.log Flag, where if do.log is true the calculations are computed on the log scale.
#'
#' @return Likelihood of the intensity function.
#' @export
#'
likelihood_x <- function(partition, heights, d.spikes, hyper.param, T.min = NULL, ISI.type = "Gamma", do.log = TRUE){

  # Firstly check everything is consistent.
  check_feasible(partition, heights, d.spikes, T.min)

  #  Get the x(y_i) and X(y_i-1,y_i) values
  values <- get_X(partition,heights,d.spikes, T.min)
  d.x <- values$intensity  ; d.X <- values$integral

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
        out <- out + sum(base::log(x[2:N]) -base::log(2*X[2:N]) - 0.5*base::log(mu*pi) - (base::log(X[2:N]) +mu)**2)/(4*mu)

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


#' likelihood of the hyper parameter.
#'
#' @param value ISI parameter.
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights  Vector of (K-1) heights corresponding to the partition.
#' @param d.spikes Data frame containing spike sequences.
#' @param prior.hyper.param Prior of the hyper parameter.
#' @param T.min Refractory period, default set to NULL.
#' @param ISI.type The ISI distribution.
#'
#' @return marginal of the hyper parameter
#' @export
log_pi_hyper_param <- function(value, partition, heights, d.spikes, prior.hyper.param, T.min = NULL, ISI.type = "Gamma"){
  # Get X for each of the spike sequences.
  d.X <- get_X(partition,heights,d.spikes,T.min)$integral

  # Initially add the value from the prior
  out <- (prior.hyper.param[1]-1) * base::log(value) - prior.hyper.param[2]*value

  ##----------------------------
  ## FOR IN GAMMA ISI DISTRIBUTION.
  ##----------------------------
  if(ISI.type == "Gamma"){
    gam <- value
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
    alpha <- value
    if (alpha < 0){
      stop("Error! Your alpha value is less then 0. alpha = ",alpha,".")
    }

    d.x <- get_X(partition,heights,d.spikes,T.min)$intensity

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
    l <- value
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
    mu <- value
    if (mu < 0){
      stop("Error! Your mu value is less then 0. mu = ",mu,".")
    }
    d.x <- get_X(partition,heights,d.spikes,T.min)$intensity

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
    k <- value
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
      }


    }
  }
  ##----------------------------
  else{
    stop("Error! You haven't chose a correct ISI.type. \n You inputted:\"", ISI.type,"\".\n")
  }

  return(out)
}

#' Likelihood of the refratory period
#'
#' @param T.min Refractory period.
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param d.spikes Data frame containing spike sequences.
#' @param hyper ISI parameter.
#' @param prior.T.min Prior for the refractory period.
#' @param ISI.type The ISI distribution.
#'
#' @return marginal of the refractory period
#' @export
#'
log_pi_Tmin <- function(T.min, partition, heights, d.spikes, hyper, prior.T.min, ISI.type = "Gamma"){
  # Get X for each of the spike sequences.
  d.X <- get_X(partition,heights,d.spikes,T.min)$integral

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
    gam <- hyper
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
    alpha <- hyper
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
    l <- hyper
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
    mu <- hyper
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
    k <- hyper
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

#' RJMCMC: Move change point.
#'
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param spikes Data frame containing spike sequences.
#' @param hyper.param ISI parameter.
#' @param T.min Refractory period.
#' @param ISI.type The ISI distribution.
#' @param do.log Flag, where if do.log is true the calculations are computed on the log scale.
#'
#' @return new partition and heights
#' @export
move_change_point <- function(partition, heights, spikes, hyper.param, T.min = NULL, ISI.type = "Gamma", do.log){
  if(length(partition) > 2){
    #  calculate the likelihood of original intensity function.
    likeli.cur <- likelihood_x(partition, heights, spikes, hyper.param, T.min, ISI.type)

    # Generate new intensity function, and calculate it's liklihood.
    M <- length(partition)
    # Use if because sample(2,size=1) samples from 1:2 rather then just 2.
    if( M==3){
      selected.pos <- 2
    }
    else{
      selected.pos <- sample(2:(M-1),size = 1)
    }
    new.position <- stats::runif(1,min = partition[selected.pos-1], max = partition[selected.pos+1])
    partition.can <- partition
    partition.can[selected.pos] <- new.position
    likeli.can <- likelihood_x(partition.can, heights, spikes, hyper.param, T.min, ISI.type)

    # Calculate the probability that he new intensity function is accepted.
    if(do.log == FALSE){
      p.acc <- (likeli.can/likeli.cur) * ( (partition[selected.pos+1] - new.position)*(new.position - partition[selected.pos-1] ) /
                                             ( (partition[selected.pos+1]-partition[selected.pos])*(partition[selected.pos]-partition[selected.pos-1])  ))

      # Draw random uniform to accept/reject new intensity function.

      u<- stats::runif(1)
      if (u < p.acc){
        partition <- partition.can
      }
    }
    else{
      p.acc <- likeli.can - likeli.cur + log( (partition[selected.pos+1] - new.position)*(new.position - partition[selected.pos-1]) ) -
        log( (partition[selected.pos+1]-partition[selected.pos])*(partition[selected.pos]-partition[selected.pos-1])  )

      u<- stats::runif(1)
      if (log(u) < p.acc){
        partition <- partition.can
      }
    }

  }
  return(partition)
}

#' RJMCMC: Change height.
#'
#' @param partition  Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param spikes Data frame containing spike sequences.
#' @param hyper.param ISI parameter.
#' @param kappa parameter of priors height
#' @param mu parameter of priors height
#' @param kappa0 parameter of priors height
#' @param T.min Refractory period.
#' @param ISI.type The ISI distribution.
#' @param do.log Flag, where if do.log is true the calculations are computed on the log scale.
#' @param which.heights Method for prior distribution of heights, either 'independent' or 'martingale'.
#'
#' @return New partition and heights
#' @export
#'
change_height <- function(partition, heights, spikes, hyper.param, kappa, mu, kappa0, T.min = NULL, ISI.type = "Gamma", do.log, which.heights = "independent"){
  #  Calculate the likelihood of original intensity function.
  likeli.cur <- likelihood_x(partition, heights, spikes, hyper.param, T.min, ISI.type, do.log)

  # Next generate new height and calculate it's likelihood.
  N <- length(heights)
  selected.pos <- base::sample(1:N, size = 1)
  v <- stats::runif(1, max = 0.5, min = -0.5)
  height.new <- heights[selected.pos] * exp(v)
  heights.can <- heights
  heights.can[selected.pos] <- height.new
  likeli.can <- likelihood_x(partition, heights.can, spikes, hyper.param, T.min, ISI.type, do.log)

  # Calculate the probability that he new intensity function is accepted.
  if(do.log == FALSE){
    # Calculate the prior on the heights.
    if(which.heights == "independent"){
      prior.heights <- (height.new/heights[selected.pos])**kappa *
        exp(-mu*(height.new - heights[selected.pos]))
    }
    else if(which.heights == "martingale"){
      if(selected.pos == N || length(heights) == 1){
        old.height <- heights[selected.pos] ;
        if(length(heights) == 1){pre.height <- mu}else{pre.height <- heights[selected.pos-1]}
        prior.heights <- (height.new/old.height)**(kappa-1)* exp(-(kappa/pre.height)*(height.new - old.height))
      }
      else{
        old.height <- heights[selected.pos] ; after.height <- heights[selected.pos+1]
        if(selected.pos == 1){pre.height <- mu}else{pre.height <- heights[selected.pos-1]}
        prior.heights <- old.height/height.new * exp(-(kappa/pre.height)*(height.new - old.height) - after.height*((kappa/height.new) - (kappa/old.height)))
      }
    }
    else{stop("Incorrect which heights!")}
    p.acc <- (likeli.can/likeli.cur) * prior.heights

    # Draw random uniform to accept/reject new intensity function.
    u<- stats::runif(1)
    if (u < p.acc){
      heights <- heights.can
    }
  }
  else{
    # Calculate the prior on the heights.
    if(which.heights == "independent"){
      prior.heights <- kappa* log((height.new/heights[selected.pos])) -
        mu*(height.new - heights[selected.pos])
    }
    else if(which.heights == "martingale"){
      if(selected.pos == N || length(heights) == 1){
        old.height <- heights[selected.pos] ;
        if(length(heights) == 1){
          prior.heights <- (kappa0-1)*(log(height.new) - log(old.height))  -(mu)*(height.new - old.height)
        }
        else{
          pre.height <- heights[selected.pos-1]
          prior.heights <- (kappa-1)*log(height.new/old.height)  -(kappa/pre.height)*(height.new - old.height)

        }
      }
      else{
        old.height <- heights[selected.pos] ; after.height <- heights[selected.pos+1]
        if(selected.pos == 1){
          prior.heights <- (kappa0-1-kappa)*(log(height.new)-log(old.height))  -(mu)*(height.new - old.height) - after.height*((kappa/height.new) - (kappa/old.height))

        }
        else{
          pre.height <- heights[selected.pos-1]
          prior.heights <- log(old.height/height.new)  -(kappa/pre.height)*(height.new - old.height) - after.height*((kappa/height.new) - (kappa/old.height))

        }
      }
    }
    else{stop("Incorrect which heights!")}

    p.acc <- (likeli.can -likeli.cur) + prior.heights

    # Draw random uniform to accept/reject new intensity function.
    u<- stats::runif(1)
    if (log(u) < p.acc){
      heights <- heights.can
    }
  }

  return(heights)
}

#' RJMCMC: Birth of a change point.
#'
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param spikes Data frame containing spike sequences.
#' @param hyper.param ISI parameter.
#' @param cvalue parameter c in the RJMCMC, to choose which event to perform in MCMC iteration.
#' @param kappa parameter of priors height
#' @param mu parameter of priors height
#' @param kappa0 parameter of priors height
#' @param k.max Maximum number of change points in RJMCMC.
#' @param sum.of.prob sum^k.max_i=0 Poi(k | lambda).
#' @param lambda Parameter of the Poisson distribution.
#' @param T.min Refractory period.
#' @param ISI.type The ISI distribution.
#' @param do.log Flag, where if do.log is true the calculations are computed on the log scale.
#' @param which.heights Method for prior distribution of heights, either 'independent' or 'martingale'.
#'
#' @return new partition and heights
#' @export
birth <- function(partition, heights, spikes, hyper.param, cvalue, kappa, mu, kappa0, k.max, sum.of.prob, lambda, T.min = NULL, ISI.type = "Gamma", do.log, which.heights = "independent"){
  #  Calculate the likelihood of original intensity function.

  likeli.cur <- likelihood_x(partition, heights, spikes, hyper.param, T.min, ISI.type, do.log)

  # Next we generate the new intensity function parameters.
  new.point <- stats::runif(1,min = partition[1], max = partition[length(partition)])
  partition.new <- sort(c(partition,new.point))
  pos.new.point <- which(abs(partition.new-new.point) < 1e-10 )
  v <- stats::runif(1, min = 0, max = 1)
  A <- exp(new.point - partition.new[pos.new.point-1] )
  B <- exp(partition.new[pos.new.point+1] -new.point)
  C <- exp(partition.new[pos.new.point+1] -partition.new[pos.new.point-1])

  new.height1 <- heights[pos.new.point-1]/ ( ((1-v)/v)**(B/C) )
  new.height2 <- ((1-v)/v) * new.height1
  if (length(partition) == 2){
    heights.new <- c(new.height1, new.height2)
  }
  else if (pos.new.point == 2){
    heights.new <- c(new.height1, new.height2, heights[2:length(heights)])
  }
  else if (pos.new.point ==length(partition.new)-1){
    heights.new <- c(heights[1:(length(heights)-1)], new.height1, new.height2)

  }
  else{
    heights.new <- c(heights[1:(pos.new.point-2)],new.height1,new.height2,heights[pos.new.point:length(heights)])
  }

  # Tie to calculate whether to accept the new function.
  # Do we do calculations on log scale?
  if(do.log == F){
    # Calculate the likelihood of new intensity function.
    likeli.can <- likelihood_x(partition.new, heights.new, spikes, hyper.param, T.min, ISI.type, do.log)

    # Next calculate the jacobian.
    jacob <- ((new.height1 + new.height2)**2 )/heights[pos.new.point-1]

    # Define k the number of change points as this is required for the proposal ratio and prior ratio.
    k <- length(partition) - 2

    # Calculate the proposal ratio, d.k1 = d_{k+1}.
    d.k1 <- cvalue* min(1, conditioned_poisson(k,lambda,k.max,sum.of.prob)/conditioned_poisson(k+1,lambda,k.max,sum.of.prob))
    b.k <- cvalue* min(1, conditioned_poisson(k+1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob))
    proposal <- ( d.k1 * (partition[length(partition)] - partition[1]) ) / (b.k * (k+1))

    # Calculate the prior ratio.
    prior.change.points <- conditioned_poisson(k+1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob)
    prior.pos.change.points <- (2*k+1)*(2*k+3)*(new.point - partition.new[pos.new.point-1] )*(partition.new[pos.new.point+1] -new.point)/
      ( ((partition[length(partition)] - partition[1])**2) * (partition.new[pos.new.point+1] -partition.new[pos.new.point-1]) )
    if(which.heights == "independent"){
      prior.heights <- (mu**kappa / gamma(kappa)) * ( (new.height1*new.height2/heights[pos.new.point-1])**(kappa-1) ) *
        (exp(-mu*(new.height1+new.height2-heights[pos.new.point-1])))
    }
    else if(which.heights == "martingale"){

      # Case new heights are at the very end or going from 1 height to 2.
      if(pos.new.point ==length(partition.new)-1 || length(heights) == 1){
        old.height <- heights[pos.new.point-1]
        if(length(heights) == 1){
          prior.heights <- (new.height1**(kappa0-1-kappa))*(kappa**kappa)*(new.height2**(kappa-1))*(old.height**(1-kappa0))*
            exp(-(mu)*(new.height1 - old.height) -(kappa/new.height1)*new.height2)
        }
        else{
          pre.height <- heights[pos.new.point-2]
          prior.heights <- ((kappa**kappa/new.height1)/gamma(kappa)) * (new.height2/old.height)**(kappa-1) *
            exp(-(kappa/pre.height)*(new.height1 - old.height) -(kappa/new.height1)*new.height2)
        }
      }
      #  Case for middle points and if new point is at the beginning.
      else{
        old.height <- heights[pos.new.point-1] ; after.height <- heights[pos.new.point]
        if(pos.new.point == 2){
          prior.heights <- (kappa**kappa / gamma(kappa)) * ((new.height1/old.height)**(kappa0-1-kappa)/new.height2)*
            exp(-(kappa/pre.height)*(new.height1 - old.height) -(kappa/new.height1)*new.height2 - after.height*((kappa/new.height2) - (kappa/old.height)))

        }
        else{
          pre.height <- heights[pos.new.point-2]

          prior.heights <- (kappa**kappa / gamma(kappa)) * (old.height/new.height1*new.height2) *
            exp(-(kappa/pre.height)*(new.height1 - old.height) -(kappa/new.height1)*new.height2 - after.height*((kappa/new.height2) - (kappa/old.height)))

        }
      }
    }
    else{
      stop("Haven't entered a correct which.heights.")
    }

    prior.ratio <- prior.change.points * prior.pos.change.points * prior.heights

    # Calculate the probability that he new intensity function is accepted.
    p.acc <- (likeli.can/likeli.cur) * prior.ratio * proposal * jacob

    # Draw random uniform to accept/reject new intensity function.
    u<- stats::runif(1)
    if (u < p.acc){
      heights <- heights.new
      partition <- partition.new
    }

  }
  else{
    # Calculations done on log scale.
    # Calculate the likelihood of new intensity function.
    likeli.can <- likelihood_x(partition.new, heights.new, spikes, hyper.param, T.min, ISI.type, do.log)

    # Next calculate the jacobian.
    jacob <- ((new.height1 + new.height2)**2 )/heights[pos.new.point-1]

    # Define k the number of change points as this is required for the proposal ratio and prior ratio.
    k <- length(partition) - 2

    # Calculate the proposal ratio, d.k1 = d_{k+1}.
    d.k1 <- cvalue* min(1, conditioned_poisson(k,lambda,k.max,sum.of.prob)/conditioned_poisson(k+1,lambda,k.max,sum.of.prob))
    b.k <- cvalue* min(1, conditioned_poisson(k+1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob))
    proposal <- ( d.k1 * (partition[length(partition)] - partition[1]) ) / (b.k * (k+1))

    # Calculate the prior ratio.
    prior.change.points <- conditioned_poisson(k+1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob)
    prior.pos.change.points <- (2*k+1)*(2*k+3)*(new.point - partition.new[pos.new.point-1] )*(partition.new[pos.new.point+1] -new.point)/
      ( ((partition[length(partition)] - partition[1])**2) * (partition.new[pos.new.point+1] -partition.new[pos.new.point-1]) )

    if(which.heights == "independent"){
      prior.heights <- kappa*log(mu) - lgamma(kappa) + (kappa-1)*(log(new.height1) + log(new.height2) - log(heights[pos.new.point-1])) -
        mu*(new.height1+new.height2-heights[pos.new.point-1])
    }
    else if(which.heights == "martingale"){

      # Case new heights are at the very end or going from 1 height to 2.
      if(pos.new.point ==length(partition.new)-1 || length(heights) == 1){
        old.height <- heights[pos.new.point-1]
        if(length(heights) == 1){
          prior.heights <- (kappa0-1-kappa)*log(new.height1) + kappa*log(kappa) + (kappa-1)*log(new.height2) - (kappa0-1)*log(old.height) -
            lgamma(kappa) - (mu)*(new.height1 - old.height) -(kappa/new.height1)*new.height2
        }
        else{
          pre.height <- heights[pos.new.point-2]
          prior.heights <- kappa*log(kappa) -lgamma(kappa) +(kappa-1)*(log(new.height2) - log(old.height)) - log(new.height1) -
            (kappa/pre.height)*(new.height1 - old.height) - (kappa/new.height1)*new.height2
        }

      }
      #  Case for middle points and if new point is at the beginning.
      else{
        old.height <- heights[pos.new.point-1] ; after.height <- heights[pos.new.point]
        if(pos.new.point == 2){
          prior.heights <- kappa*log(kappa) - lgamma(kappa) + (kappa0-1-kappa)*(log(new.height1)-log(old.height)) - log(new.height2)-
            mu*(new.height1- old.height) - kappa* ((after.height/new.height2) - (after.height/old.height) + (new.height2/new.height1))
        }
        else{
          pre.height <- heights[pos.new.point-2]
          prior.heights <- kappa*log(kappa) - lgamma(kappa) +log(old.height)-log(new.height2)-log(new.height1) -
            (kappa/pre.height)*(new.height1 - old.height) - kappa*((after.height/new.height2) - (after.height/old.height) + (new.height2/new.height1))

        }
      }
    }
    else{
      stop("Haven't entered a correct which.heights.")
    }

    prior.ratio <- log(prior.change.points) + log(prior.pos.change.points) + prior.heights

    p.acc <- (likeli.can-likeli.cur) + prior.ratio + log(proposal) + log(jacob)

    # Draw random uniform to accept/reject new intensity function.
    u<- stats::runif(1)
    if (log(u) < p.acc){
      heights <- heights.new
      partition <- partition.new
    }


  }

  # return(p.acc)
  return(list(partition = partition, heights = heights))
}

#' RJMCMC: Death of a change point
#'
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param spikes Data frame containing spike sequences.
#' @param hyper.param ISI parameter.
#' @param cvalue parameter c in the RJMCMC, to choose which event to perform in MCMC iteration.
#' @param kappa parameter of priors height
#' @param mu parameter of priors height
#' @param kappa0 parameter of priors height
#' @param k.max Maximum number of change points in RJMCMC.
#' @param sum.of.prob sum^k.max_i=0 Poi(k | lambda).
#' @param lambda Parameter of the Poisson distribution.
#' @param T.min Refractory period.
#' @param ISI.type The ISI distribution.
#' @param do.log Flag, where if do.log is true the calculations are computed on the log scale.
#' @param which.heights Method for prior distribution of heights, either 'independent' or 'martingale'.
#'
#' @return new partition and heights
#' @export
#'
#' @examples
death <- function(partition, heights, spikes, hyper.param, cvalue, kappa, mu, kappa0, k.max, sum.of.prob, lambda, T.min = NULL, ISI.type = "Gamma", do.log, which.heights = "independent"){
  N <- length(partition)

  # If to make sure there exists a point to delete.
  if (N >2){
    #  Calculate the likelihood of original intensity function.
    likeli.cur <- likelihood_x(partition, heights, spikes, hyper.param, T.min, ISI.type, do.log)

    # Next we generate the new intensity function parameters.
    if (N == 3){pos.to.die <- 2}
    else{
      pos.to.die <- base::sample(2:(N-1),1)
    }
    A <- partition[pos.to.die] - partition[pos.to.die-1]
    B <- partition[pos.to.die+1] - partition[pos.to.die]
    C <- partition[pos.to.die+1] - partition[pos.to.die-1]
    new.height <- exp( (A*log(heights[pos.to.die-1]) + B*log(heights[pos.to.die]))/C )

    partition.new <- partition[-pos.to.die]
    if (length(partition.new) == 2){
      heights.new <- new.height
    }
    else if (pos.to.die == 2){
      heights.new <- c(new.height, heights[3:(N-1)])
    }
    else if (pos.to.die == N-1){
      heights.new <- c(heights[1:(N-3)], new.height)
    }
    else{
      heights.new <- c(heights[1:(pos.to.die-2)], new.height, heights[(pos.to.die+1):length(heights)] )
    }

    #Do we calculate on the log scale.
    if(do.log== F){
      # Calculate the likelihood of new intensity function.
      likeli.can <- likelihood_x(partition.new, heights.new, spikes, hyper.param, T.min, ISI.type, do.log)

      #  --------------------------------------------------
      # Next calculate the jacobian.
      jacob <- new.height / ( (heights[pos.to.die] + heights[pos.to.die-1])**2 )

      # Define k the number of change points as this is required for the proposal ratio and prior ratio.
      k <- length(partition) - 2

      # Calculate the proposal ratio, b.k1 = b_{k-1}.
      d.k <- cvalue* min(1, conditioned_poisson(k-1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob))
      b.k1 <- cvalue* min(1, conditioned_poisson(k,lambda,k.max,sum.of.prob)/conditioned_poisson(k-1,lambda,k.max,sum.of.prob))
      proposal <- (b.k1 * k) / ( d.k * (partition[length(partition)] - partition[1]) )

      # Calculate the prior ratio.
      prior.change.points <- conditioned_poisson(k-1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob)
      prior.pos.change.points <- (2*k+1)*(2*k)*(partition[pos.to.die+1] - partition[pos.to.die-1] )/
        ( ((partition[length(partition)] - partition[1])**2) * (partition[pos.to.die] -partition[pos.to.die-1]) * (partition[pos.to.die+1]-partition[pos.to.die]) )
      if(which.heights == "independent"){
        prior.heights <- (gamma(kappa)/ mu**kappa) * ( (new.height/(heights[pos.to.die]*heights[pos.to.die-1]) )**(kappa-1) ) *
          (exp(-mu*(new.height-heights[pos.to.die-1]-heights[pos.to.die])))
      }
      else if(which.heights == "martingale"){

        # Case new heights are at the very end or going from 1 height to 2.
        if(pos.to.die ==length(heights) || length(heights) == 2){
          old.height1 <- heights[pos.to.die-1] ;  old.height2 <- heights[pos.to.die] ; after.height <- heights[pos.to.die+1]
          if(length(heights) == 2){
            prior.heights <- gamma(kappa)/( (kappa**kappa)*(old.height1**(kappa0-1-kappa))*(old.height2**(kappa-1))*(new.height**(kappa0-1)) ) +
              mu*(old.height1-new.height) + kappa*old.height2/old.height1
          }
          else{
            pre.height <- heights[pos.to.die-2]
            prior.heights <- (old.height1*gamma(kappa)/ kappa**kappa) * (new.height/old.height2)**(kappa-1) *
              exp((kappa/pre.height)*(old.height1 - new.height) +(kappa/old.height1)*old.height2)
          }

        }
        #  Case for middle points and if new point is at the beginning.
        else{
          old.height1 <- heights[pos.to.die-1] ;  old.height2 <- heights[pos.to.die] ; after.height <- heights[pos.to.die+1]
          if(pos.to.die == 2){
            prior.heights <- gamma(kappa) *((new.height/old.height1)**(kappa0-1-kappa))*after.height *
              exp(mu*(new.height-old.height1) - kappa*(after.height/new.height + old.height2/old.height1 + after.height/old.height1 )) / (kappa**kappa)
          }
          else{
            pre.height <- heights[pos.to.die-2]
            prior.heights <- (gamma(kappa) / kappa**kappa) * (old.height1*old.height2/new.height) *
              exp((kappa/pre.height)*(old.height1 - new.height) +(kappa/old.height1)*old.height2 + after.height*((kappa/old.height2) - (kappa/new.height)))

          }
        }
      }
      else{
        stop("Haven't entered a correct which.heights.")
      }
      prior.ratio <- prior.change.points * prior.pos.change.points * prior.heights

      #  --------------------------------------------------

      # Calculate the probability that he new intensity function is accepted.
      p.acc <- (likeli.can/likeli.cur) * prior.ratio * proposal * jacob

      # Draw random uniform to accept/reject new intensity function.
      u<- stats::runif(1)
      if (u < p.acc){
        heights <- heights.new
        partition <- partition.new
      }
    }
    else{
      # Calculate the likelihood of new intensity function.
      likeli.can <- likelihood_x(partition.new, heights.new, spikes, hyper.param, T.min, ISI.type, do.log)

      #  --------------------------------------------------
      # Next calculate the jacobian.
      jacob <- new.height / ( (heights[pos.to.die] + heights[pos.to.die-1])**2 )

      # Define k the number of change points as this is required for the proposal ratio and prior ratio.
      k <- length(partition) - 2

      # Calculate the proposal ratio, b.k1 = b_{k-1}.
      d.k <- cvalue* min(1, conditioned_poisson(k-1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob))
      b.k1 <- cvalue* min(1, conditioned_poisson(k,lambda,k.max,sum.of.prob)/conditioned_poisson(k-1,lambda,k.max,sum.of.prob))
      proposal <- (b.k1 * k) / ( d.k * (partition[length(partition)] - partition[1]) )

      # Calculate the prior ratio.
      prior.change.points <- conditioned_poisson(k-1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob)
      prior.pos.change.points <- (2*k+1)*(2*k)*(partition[pos.to.die+1] - partition[pos.to.die-1] )/
        ( ((partition[length(partition)] - partition[1])**2) * (partition[pos.to.die] -partition[pos.to.die-1]) * (partition[pos.to.die+1]-partition[pos.to.die]) )
      if(which.heights == "independent"){
        prior.heights <- lgamma(kappa) - kappa*log(mu) + (kappa-1)*(log(new.height) - log(heights[pos.to.die]) - log(heights[pos.to.die-1])) -
          mu*(new.height-heights[pos.to.die-1]-heights[pos.to.die])
      }
      else if(which.heights == "martingale"){

        # Case new heights are at the very end or going from 1 height to 2.
        if(pos.to.die ==length(heights) || length(heights) == 2){
          old.height1 <- heights[pos.to.die-1] ;  old.height2 <- heights[pos.to.die] ; after.height <- heights[pos.to.die+1]
          if(length(heights) == 2){
            prior.heights <- lgamma(kappa) - kappa*log(kappa) - (kappa0-1-kappa)*log(old.height1) - (kappa-1)*log(old.height2) +
              (kappa0-1)*log(new.height) -mu*(new.height-old.height1) + kappa*old.height2/old.height1
          }
          else{
            pre.height <- heights[pos.to.die-2]
            prior.heights <- lgamma(kappa) +log(old.height1) - kappa*log(kappa) +(kappa-1)*(log(new.height)-log(old.height2)) +
              (kappa/pre.height)*(old.height1 - new.height) +(kappa/old.height1)*old.height2
          }
        }
        #  Case for middle points and if new point is at the beginning.
        else{
          old.height1 <- heights[pos.to.die-1] ;  old.height2 <- heights[pos.to.die] ; after.height <- heights[pos.to.die+1]
          if(pos.to.die == 2){
            prior.heights <- lgamma(kappa) + (kappa0-1-kappa)*(log(new.height)-log(old.height1)) +log(old.height2) - kappa*log(kappa) -
              mu*(new.height-old.height1) -kappa*( after.height/new.height - after.height/old.height2 - old.height2/old.height1)
          }
          else{
            pre.height <- heights[pos.to.die-2]
            prior.heights <- lgamma(kappa) - kappa*log(kappa) +log(old.height1) +log(old.height2) - log(new.height) -
              (kappa/pre.height)*(new.height - old.height1) - kappa*(after.height/new.height - old.height2/old.height1 - after.height/old.height2)

          }
        }
      }
      else{
        stop("Haven't entered a correct which.heights.")
      }
      prior.ratio <- log(prior.change.points) + log(prior.pos.change.points) + prior.heights

      #  --------------------------------------------------

      # Calculate the probability that he new intensity function is accepted.
      p.acc <- (likeli.can-likeli.cur) + prior.ratio + log(proposal) + log(jacob)

      # Draw random uniform to accept/reject new intensity function.
      u<- stats::runif(1)
      if (log(u) < p.acc){
        heights <- heights.new
        partition <- partition.new
      }

    }
  }
  # return(p.acc)
  return(list(partition = partition, heights = heights))
}


#' RJMCMC: MCMC.
#' @param spikes Data frame containing spike sequences.
#' @param end.time End time of the experiment
#' @param iter Number of iterations to record
#' @param burn Number of iterations to burn before recording.
#' @param k.max Maximum number of change points in RJMCMC.
#' @param lambda Parameter of the Poisson distribution.
#' @param kappa parameter of priors height
#' @param mu parameter of priors height
#' @param kappa0 parameter of priors height
#' @param hyper.param Value of the hyper parameter, or values for the gamma prior of the hyper parameter.
#' @param sigma.h Sigma used in RW-Metropolis for the ISI parameter.
#' @param start.hyper Iteration number to begin the ISI parameter.
#' @param hyper.initial Initial value of the ISI parameter.
#' @param T.min.param Value of the refractory period, or values for the gamma prior of the refractory period parameter.
#' @param T.min.initial Initial value of the refractory period.
#' @param sigma.t Sigma used in RW-Metropolis for the refractory period parameter.
#' @param ISI.type The ISI distribution.
#' @param do.log Flag, where if do.log is true the calculations are computed on the log scale.
#' @param show.iter Flag for whether to print the number of iterations complete into the console.
#' @param which.heights Method for prior distribution of heights, either 'independent' or 'martingale'.
#'
#' @return Iterations of the MCMC algorithm.
#' @export
#'
mcmc <- function(spikes,end.time,iter,burn, k.max, lambda, kappa, mu, kappa0,
                 hyper.param = c(1,0.001), sigma.h = NULL, start.hyper = 1000, hyper.initial = 0.1,
                 T.min.param = NULL, T.min.initial = 0.05, sigma.t = NULL, ISI.type = "Gamma",
                 do.log = TRUE, show.iter = FALSE, which.heights = "independent"){

  # Print inputs if to check.
  # print("HI we in mcmc")
  #  print(spikes)
  # print(end.time)
  # print(iter)
  # print(burn)
  # print(k.max)
  # print(lambda)
  # print(kappa)
  # print(mu)


  # Check gammma.param has length 1 or 2.
  if(length(hyper.param)>2){stop("Error! The length of hyper.param != 1 or 2.\n The hyper.param you entered has length ",
                                 length(hyper.param),". \n")}
  # Check T.min.param has length 1 or 2.
  if(length(T.min.param)>2){stop("Error! The length of T.min.param != 1 or 2.\n The hyper.param you entered has length ",
                                 length(T.min.param),". \n")}

  # Create output array to store the iterations, for x(t), and vectors for hyper and T.min if required.
  out.part <- matrix(NA, nrow = iter, ncol = k.max+2)
  out.heights <- matrix(NA,nrow = iter, ncol = k.max+1)
  if(length(hyper.param)==2){
    out.hyper <- rep(NA,iter)
    if(is.null(sigma.h) == TRUE){stop("Error! If you want to use mcmc to get ISI parameter estimation you require a sigma.h value.")}
  }
  if(length(T.min.param)==2){
    out.T.min <- rep(NA,iter)
    if(is.null(sigma.t) == TRUE){stop("Error! If you want to use mcmc to get T.min estimation you require a sigma.t value.")}

    # Define T.max so that we don't put forward a value larger then possible.
    T.max <- 100000
    for(i in 1:ncol(spikes)){
      s <- spikes[,i]
      s <- s[!is.na(s)]
      new.T.max <- min(s[-1] - s[-length(s)])
      T.max <- min(T.max,new.T.max)
    }
  }

  # Calculate the sum of prob, to be used in calculating the probability of k steps.
  sum.of.prob <- sum(stats::dpois(0:k.max,lambda))

  # Next calcuate the value of c.
  val <- -10
  for ( i in 0:k.max){
    val <- max(val, min(1, conditioned_poisson(i+1,lambda,k.max,sum.of.prob)/conditioned_poisson(i,lambda,k.max,sum.of.prob)) +
                 min(1,conditioned_poisson(i,lambda,k.max,sum.of.prob)/conditioned_poisson(i-1,lambda,k.max,sum.of.prob) ))
  }
  cvalue <- 0.9 / val

  # setup the initial step function and gamma values.
  partition <- c(0,end.time)
  heights <- length(spikes)/end.time
  if (length(hyper.param)==2){hyper.cur <- hyper.initial}
  else {hyper.cur <- hyper.param}
  if (length(T.min.param)==2){T.min.cur <- T.min.initial}
  else {T.min.cur <- T.min.param}

  # mcmc loop starts here
  li<- seq(100,150000,100)
  for (i in 2:(iter+burn)) {
    # Print the loop.
    if (i %in% li && show.iter == TRUE){ cat(i,"\n")}

    #############################
    # UPDATE INTENSITY FUNCTION
    #############################
    # Draw a uniform r.v to calculate what to update
    u <- stats::runif(1)
    k <- length(partition) - 2
    b.k <- cvalue * min(1, conditioned_poisson(k+1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob))
    d.k <- cvalue * min(1, conditioned_poisson(k-1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob))


    if (u < b.k){
      next.val <- birth(partition, heights, spikes, hyper.cur, cvalue, kappa, mu, kappa0, k.max, sum.of.prob, lambda, T.min.cur, ISI.type, do.log, which.heights)
      partition <- next.val$partition
      heights <- next.val$heights
    }
    else if (b.k < u  && u < b.k + d.k){
      next.val <- death(partition, heights, spikes, hyper.cur, cvalue, kappa, mu, kappa0, k.max, sum.of.prob, lambda, T.min.cur, ISI.type, do.log, which.heights)
      partition <- next.val$partition
      heights <- next.val$heights
    }
    else if (b.k + d.k < u && u < 1){
      partition <- move_change_point(partition, heights, spikes, hyper.cur, T.min.cur, ISI.type, do.log)
      heights <- change_height(partition, heights, spikes, hyper.cur, kappa, kappa0, mu, T.min.cur, ISI.type, do.log, which.heights)
    }
    else{ stop("Error A1! Something gone wrong!")}


    #############################
    # UPDATE HYPER PARAMETERS
    #############################
    # if(i==100){
    #   ratio <- length(spikes)/end.time
    #   print(ratio)
    # }
    if(length(hyper.param)==2 && i > start.hyper){
      # print(hyper.cur)
      hyper.can <- stats::rnorm(1, hyper.cur, sigma.h)
      hyper.can
      if (hyper.can > 0) {
        # print(hyper.cur)
        # print( partition)
        # print(spikes)
        # print( hyper.param)
        # print(T.min.cur)
        # print(ISI.type)
        # print(heights)
        log.pi.cur <- log_pi_hyper_param(hyper.cur, partition, heights, spikes, hyper.param, T.min.cur, ISI.type)
        log.pi.can <- log_pi_hyper_param(hyper.can, partition, heights, spikes, hyper.param, T.min.cur, ISI.type)
        # M-H ratio
        # cat("cur = ", log.pi.cur,"\n")
        # cat("can = ", log.pi.can,"\n")
        # draw from a U(0,1)
        u <- stats::runif(1)

        if (log(u) < log.pi.can - log.pi.cur) {
          hyper.cur <- hyper.can
        }
      }
    }

    #############################
    # UPDATE T.min
    #############################

    if(length(T.min.param)==2){
      T.min.can <- stats::rnorm(1, T.min.cur, sigma.t)
      if (T.min.can > 0 && T.min.can < T.max) {
        log.pi.cur <- log_pi_Tmin(T.min.cur, partition, heights, spikes, hyper.cur, T.min.param, ISI.type)
        log.pi.can <- log_pi_Tmin(T.min.can, partition, heights, spikes, hyper.cur, T.min.param, ISI.type)
        # M-H ratio

        # draw from a U(0,1)
        u <- stats::runif(1)

        if (log(u) < log.pi.can - log.pi.cur) {
          T.min.cur <- T.min.can
        }
      }
    }



    #############################
    # STORE THE OUTPUTS
    #############################
    if(i > burn){
      j<- i-burn
      out.part[j,1:length(partition)] <- partition
      out.heights[j,1:length(heights)] <- heights
      if(length(hyper.param)==2){ out.hyper[j] <- hyper.cur}
      if(length(T.min.param)==2){ out.T.min[j] <- T.min.cur}
    }

  }
  hyp <- NULL ; t <- NULL
  intensity <- list(partition = out.part, heights = out.heights)
  if (length(hyper.param)==2){ hyp <- list(hyper = out.hyper)}
  if (length(T.min.param)==2){ t <- list(T.min = out.T.min)}
  # print(intensity)
  # print(hyp)
  # print(t)
  output <- do.call("c",list(intensity,hyp,t))
  return(output)
}

