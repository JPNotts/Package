#' Poisson conditioned on max value.
#'
#' @param k k
#' @param lambda Parameter of the Poisson distribution.
#' @param k.max Maximum number of change points in RJMCMC.
#' @param sum.of.prob sum^k.max_i=0 Poi(k | lambda).
#'
#' @return probability of k
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


#' RJMCMC: Move change point.
#'
#' @param partition Partition of time experiment time. A vector of K values where the first entry is 0 and the last is the end of experiment time.
#' @param heights Vector of (K-1) heights corresponding to the partition.
#' @param d.spikes Data frame containing spike sequences.
#' @param hyper.param ISI parameter.
#' @param T.min Refractory period.
#' @param ISI.type The ISI distribution.
#' @param do.log Flag, where if do.log is true the calculations are computed on the log scale.
#' @param max.T.min Maximal allowed T.min values from d.spikes.
#'
#' @return new partition and heights
move_change_point <- function(partition, heights, d.spikes, hyper.param, T.min = NULL, max.T.min = NA, ISI.type = "Gamma", do.log){
  if(length(partition) > 2){
    #  calculate the likelihood of original intensity function.
    likeli.cur <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition, heights = heights, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
    likeli.can <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition.can, heights = heights, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
#' @inheritParams move_change_point
#' @param kappa parameter of priors height
#' @param mu parameter of priors height
#' @param kappa0 parameter of priors height
#' @param which.heights Method for prior distribution of heights, either 'independent' or 'martingale'.
#'
#' @return New partition and heights
#'
change_height <- function(partition, heights, d.spikes, hyper.param, kappa, mu, kappa0, T.min = NULL, max.T.min = NA, ISI.type = "Gamma", do.log, which.heights = "independent"){
  #  Calculate the likelihood of original intensity function.
  likeli.cur <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition, heights = heights, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

  # Next generate new height and calculate it's likelihood.
  N <- length(heights)
  selected.pos <- base::sample(1:N, size = 1)
  v <- stats::runif(1, max = 0.5, min = -0.5)
  height.new <- heights[selected.pos] * exp(v)
  heights.can <- heights
  heights.can[selected.pos] <- height.new
  likeli.can <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition, heights = heights.can, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
#' @inheritParams change_height
#' @inheritParams conditioned_poisson
#' @param cvalue parameter c in the RJMCMC, to choose which event to perform in MCMC iteration.
#'
#' @return new partition and heights
birth <- function(partition, heights, d.spikes, hyper.param, cvalue, kappa, mu, kappa0, k.max, sum.of.prob, lambda, T.min = NULL, max.T.min = NA, ISI.type = "Gamma", do.log, which.heights = "independent"){
  #  Calculate the likelihood of original intensity function.

  likeli.cur <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition, heights = heights, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
    likeli.can <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition.new, heights = heights.new, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
    likeli.can <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition.new, heights = heights.new, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
#' @inheritParams birth
#'
#' @return new partition and heights
#'
#' @examples
death <- function(partition, heights, d.spikes, hyper.param, cvalue, kappa, mu, kappa0, k.max, sum.of.prob, lambda, T.min = NULL, max.T.min = NA, ISI.type = "Gamma", do.log, which.heights = "independent"){
  N <- length(partition)

  # If to make sure there exists a point to delete.
  if (N >2){
    #  Calculate the likelihood of original intensity function.
    likeli.cur <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition, heights = heights, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
      likeli.can <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition.new, heights = heights.new, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
      likeli.can <- log_pi_x(d.spikes = d.spikes, hyper.param = hyper.param, partition = partition.new, heights = heights.new, T.min = T.min, max.T.min = max.T.min,  ISI.type = ISI.type, do.log = do.log)

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
#' @param d.spikes Data frame containing spike sequences.
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
mcmc_pwc <- function(d.spikes,end.time,iter,burn, k.max, lambda, kappa, mu, kappa0,
                 hyper.param = c(1,0.001), sigma.h = NULL, start.hyper = 1000, hyper.initial = 0.1,
                 T.min.param = NULL, T.min.initial = 0.05, sigma.t = NULL, ISI.type = "Gamma",
                 do.log = TRUE, show.iter = FALSE, which.heights = "independent"){

  # Print inputs if to check.
  # print("HI we in mcmc")
  #  print(d.spikes)
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
  out.likeli <- rep(NA,iter)
  if(length(hyper.param)==2){
    out.hyper <- rep(NA,iter)
    if(is.null(sigma.h) == TRUE){stop("Error! If you want to use mcmc to get ISI parameter estimation you require a sigma.h value.")}
  }
  if(length(T.min.param)==2){
    out.T.min <- rep(NA,iter)
    if(is.null(sigma.t) == TRUE){stop("Error! If you want to use mcmc to get T.min estimation you require a sigma.t value.")}

    # Define T.max so that we don't put forward a value larger then possible.
    max.T.min <- 100000
    for(i in 1:ncol(d.spikes)){
      s <- d.spikes[,i]
      s <- s[!is.na(s)]
      new.T.max <- min(s[-1] - s[-length(s)])
      max.T.min <- min(max.T.min,new.T.max)
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
  heights <- length(d.spikes)/end.time
  if (length(hyper.param)==2){hyper.cur <- hyper.initial}
  else {hyper.cur <- hyper.param}
  if (length(T.min.param)==2){T.min.cur <- T.min.initial}
  else {T.min.cur <- T.min.param}

  # mcmc loop starts here
  pb <- utils::txtProgressBar(min = 1, max = iter+burn, style = 3)
  for (i in 2:(iter+burn)) {
    # Print the loop.
    # if (i %in% li && show.iter == TRUE){ cat(i,"\n")}

    #############################
    # UPDATE INTENSITY FUNCTION
    #############################
    # Draw a uniform r.v to calculate what to update
    u <- stats::runif(1)
    k <- length(partition) - 2
    b.k <- cvalue * min(1, conditioned_poisson(k+1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob))
    d.k <- cvalue * min(1, conditioned_poisson(k-1,lambda,k.max,sum.of.prob)/conditioned_poisson(k,lambda,k.max,sum.of.prob))


    if (u < b.k){
      next.val <- birth(partition=partition, heights=heights, d.spikes=d.spikes, hyper.param=hyper.cur, cvalue=cvalue, kappa=kappa, mu=mu, kappa0=kappa0, k.max=k.max, sum.of.prob=sum.of.prob, lambda=lambda, T.min=T.min.cur,max.T.min=max.T.min, ISI.type=ISI.type, do.log=do.log, which.heights=which.heights)
      partition <- next.val$partition
      heights <- next.val$heights
    }
    else if (b.k < u  && u < b.k + d.k){
      next.val <- death(partition=partition, heights=heights, d.spikes=d.spikes, hyper.param=hyper.cur, cvalue =cvalue, kappa=kappa, mu=mu, kappa0=kappa0, k.max=k.max, sum.of.prob=sum.of.prob, lambda=lambda, T.min=T.min.cur, max.T.min=max.T.min, ISI.type=ISI.type, do.log=do.log, which.heights=which.heights)
      partition <- next.val$partition
      heights <- next.val$heights
    }
    else if (b.k + d.k < u && u < 1){
      partition <- move_change_point(partition=partition, heights=heights, d.spikes=d.spikes, hyper.param=hyper.cur, T.min=T.min.cur, max.T.min = max.T.min, ISI.type=ISI.type, do.log=do.log)
      heights <- change_height(partition=partition, heights=heights, d.spikes=d.spikes, hyper.param=hyper.cur, kappa=kappa, kappa0=kappa0, mu=mu, T.min=T.min.cur, max.T.min=max.T.min, ISI.type=ISI.type, do.log=do.log, which.heights=which.heights)
    }
    else{ stop("Error A1! Something gone wrong!")}


    #############################
    # UPDATE HYPER PARAMETERS
    #############################
    # if(i==100){
    #   ratio <- length(d.spikes)/end.time
    #   print(ratio)
    # }
    if(length(hyper.param)==2 && i > start.hyper){
      # print(hyper.cur)
      hyper.can <- stats::rnorm(1, hyper.cur, sigma.h)
      hyper.can
      if (hyper.can > 0) {
        # print(hyper.cur)
        # print( partition)
        # print(d.spikes)
        # print( hyper.param)
        # print(T.min.cur)
        # print(ISI.type)
        # print(heights)
        log.pi.cur <- log_pi_hyper_param(d.spikes = d.spikes, hyper.param = hyper.cur, partition = partition, heights = heights,  prior.hyper.param = hyper.param, T.min = T.min.cur, max.T.min = max.T.min, ISI.type = ISI.type)
        log.pi.can <- log_pi_hyper_param(d.spikes = d.spikes, hyper.param = hyper.can, partition = partition, heights = heights,  prior.hyper.param = hyper.param, T.min = T.min.cur, max.T.min = max.T.min, ISI.type = ISI.type)
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
      if (T.min.can > 0 && T.min.can < max.T.min) {




        log.pi.cur <- log_pi_Tmin(d.spikes = d.spikes, T.min = T.min.cur, max.T.min = max.T.min, hyper = hyper.cur, partition = partition, heights = heights, prior.T.min = T.min.param, ISI.type =ISI.type)
        log.pi.can <- log_pi_Tmin(d.spikes = d.spikes, T.min = T.min.can, max.T.min = max.T.min, hyper = hyper.cur, partition = partition, heights = heights, prior.T.min = T.min.param, ISI.type =ISI.type)
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
      out.likeli[j] <- log_pi_x(d.spikes, hyper.cur, partition = partition, heights= heights, T.min = T.min.cur, max.T.min = max.T.min, ISI.type = ISI.type, do.log = TRUE)
    }
    utils::setTxtProgressBar(pb, i) # update text progress bar after each iter

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

#' get the mean and 95% conf interval for PWC method
#'
#' @param partitions  za
#' @param heights za
#' @param x za
#' @param debug za
#' @param return.val za
#' @param method za
#' @param steps za
#' @param incProg za
#'
#' @return A list
#' @export
#'
#' @examples
conf_interval_pwc <- function(partitions, heights, x, debug = FALSE, return.val = FALSE, method = "TRUE", steps = 1000, incProg = FALSE){

  # First if all iterations are constant functions
  if(dim(partitions)[2] == 2){
    quant <- quantile(heights,probs <- c(0.025,0.975))
    return(list(mean = mean(heights), lower = quant[1], upper = quant[2], partition = partitions[1,]))
  }

  # Next we consider the case where we want the exact answer (can take a while to run)
  if(method == "TRUE"){
    iter <- dim(partitions)[1]
    ex <- partitions[1,]
    ex <- ex[!is.na(ex)]
    end.time <- ex[length(ex)]
    # Create a vector of all the step postitions.
    all.step <- sort(unique(partitions[!is.na(partitions)] ))

    # Create a vector to store the mean and matrix to store the confidence interval values.
    mean.vect <- rep(NULL, length(all.step)-1)
    conf.inter <- matrix(NA,nrow =2, ncol = length(all.step)-1)

    # Compute mean and confidence interval values for the first step.
    heights.in.step <- heights[,1]
    quant <- quantile(heights.in.step,probs <- c(0.025,0.975))
    mean.vect[1] <- mean(heights.in.step)
    conf.inter[1,1] <- quant[1]
    conf.inter[2,1] <- quant[2]

    # Compute mean and confidence interval values for all other steps.
    for (i in 2:(length(all.step)-1)){
      heights.in.step <- NULL
      for(j in 1:iter){
        position <-  max(which(partitions[j,]< (all.step[i] +1e-10)))
        heights.in.step <- c(heights.in.step,heights[j,position])
      }
      quant <- quantile(heights.in.step,probs <- c(0.025,0.975))
      mean.vect[i] <- mean(heights.in.step)
      conf.inter[1,i] <- quant[1]
      conf.inter[2,i] <- quant[2]
    }

    if(return.val == TRUE){
      return(list(mean = mean.vect, lower = conf.inter[1,], upper = conf.inter[2,], partition = all.step))
    }
    else{
      # Create the plot.
      plot(0,0,ylim=c(0,(max(heights[!is.na(heights)])+0.1)), xlim = c(0,end.time), type="n", ylab = "intensity", xlab="time")

      # plot the mean value. (step function separate)
      # for( j in 1:length(mean.vect)){
      #   height <- mean.vect[j]
      #   corresp.points <- c(all.step[j],all.step[j+1])
      #   lines(corresp.points,c(height,height), col="red", lwd=2)
      # }

      # Plot mean value (step function joined)
      height.grid <- NULL
      corresp.points <- NULL
      for( j in 1:length(mean.vect)){
        height <- mean.vect[j]
        corresp.points <- c(corresp.points, all.step[j],all.step[j+1])
        height.grid <- c(height.grid, height,height)

      }
      lines(corresp.points,height.grid, col="black", lwd=2)

      # Plot the confidence interval.
      lv.grid <- NULL
      uv.grid <- NULL
      corresp.points <- NULL
      for( j in 1:length(mean.vect)){
        lv <- conf.inter[1,j]
        uv <- conf.inter[2,j]
        corresp.points <- c(corresp.points, all.step[j],all.step[j+1])
        lv.grid <- c(lv.grid, lv,lv)
        uv.grid <- c(uv.grid, uv,uv)
      }
      polygon(c(corresp.points,rev(corresp.points)), c(uv.grid, rev(lv.grid)),col = rgb(85,75,85,max=255,alpha=100),border=NA)

      # Plot the truth (if we know it)
      if (!is.null(x)){
        grid <- seq(0,end.time,end.time/(length(x)-1))
        lines(grid,x,col="red", lwd=3)
      }

    }

    if (debug == TRUE){
      cat("partitions = \n")
      print(partitions)
      cat("all.step = \n")
      print(all.step)
      cat("\n")
      cat("heights = \n")
      print(heights)
      cat("mean.vect = \n")
      print(mean.vect)
      cat("conf.inter = \n")
      print(conf.inter)
    }
  }

  # Then consider a quicker version which calculates an approximate answer.
  if(method == "Approx"){
    iter <- dim(partitions)[1]
    ex <- partitions[1,]
    ex <- ex[!is.na(ex)]
    end.time <- ex[length(ex)]
    t <- seq(0,end.time,end.time/steps)

    # Create a store for each of the outputs
    mean <- rep(NA,steps+1)
    lower <- rep(NA,steps+1)
    upper <- rep(NA,steps+1)
    # Loop over all the time points where we want to know the intensity
    store <- rep(NA,iter)
    current <- rep(1,iter)
    for(i in 1:(length(t)-1)){


      for(j in 1:iter){
        end <- sum(!is.na(partitions[j,]))
        current[j] <-  max(which(partitions[j,1:end]< (t[i] +1e-10) ))
        store[j] <- heights[j,current[j]]
      }
      mean[i] <- mean(store)
      quant <- quantile(store,probs <- c(0.025,0.975))
      lower[i] <- quant[1] ; upper[i] <- quant[2]
    }
    mean[steps+1] <- mean(store); lower[steps+1] <- quant[1] ; upper[steps+1] <- quant[2]

    if(return.val){
      return(list(mean = mean, lower = lower, upper = upper))
    }
    else{
      part <- partitions[1,]
      part <- part[!is.na(part)]
      end.time <- part[length(part)]

      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 1, 0, 0))
      plot(0,0, type ="n", xlim = c(0,end.time),
           ylim = c(min(lower) , max(upper)),
           xlab ="Time ($s$)", ylab = 'Intensity (spikes/s)', cex.lab = 1.5, cex.axis =1.2)

      t <- seq(0,end.time,end.time/1000)
      polygon(c(t,rev(t)), c(upper, rev(lower)),col = rgb(0,0,255,max=255,alpha=80),border=NA)
      lines(t,mean, col= "blue", lwd = 2)

      if(!is.null(x)){
        t<- seq(0,end.time,end.time/(length(x)-1))
        lines(t,x,col = 1, lwd = 2)
      }
    }


  }

}
