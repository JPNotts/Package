# GP inference


#' Square Exponential covariance function
#'
#' @param t  Time difference
#' @param sigma Signal Variance
#' @param l Length Scale
#'
#' @return The Covariance
#' @export
#'
cov.fn <- function(t, sigma = sqrt(1000), l = 3) {
  sigma^2 * exp(-t^2 / (2 * l^2))
}

#' Simpson's rule
#'
#' @param x x
#' @param y y
#'
#' @return integral?
#' @export
#'
cumsimpson <- function(x, y) {
  ry <- 1L
  cy <- length(y)

  dx <- diff(x)
  dy <- diff(y)
  dx1 <- dx[-length(dx)]
  dx2 <- dx[-1L]
  dxs <- dx1 + dx2
  dy1 <- dy[-length(dy)]
  dy2 <- dy[-1L]
  a <- (dy2 / (dx2 * dxs) - dy1 / (dx2 * dxs)) / 3
  b <- (dy2 * dx1 / (dx2 * dxs) + dy1 * dx2 / (dx1 * dxs)) / 2
  c <- y[-c(1L, length(y))]

  i1 <- ((a * dx1 - b) * dx1 + c) * dx1 # left half integral
  i2 <- ((a * dx2 + b) * dx2 + c) * dx2 # right half integral

  z <- vector("numeric", length(y))
  z[seq(2, length(z) - 1, 2)] <- i1[seq(1, length(i1), 2)]
  z[seq(3, length(z), 2)] <- i2[seq(1, length(i2), 2)]
  z[length(z)] <- i2[length(i2)]
  cumsum(z)
}

#' Simulate a Gaussian Process using Spectral decomposition
#'
#' @param cov.fn  The covariance function
#' @param l the length scale
#' @param dt.GP The step size
#' @param steps Number of steps required
#'
#' @return A draw from the Gaussian process
#' @export
#'
#' @examples
stat_gp <- function(cov.fn, l,dt.GP = 0.8, steps = 1000) {
  ## Generates a GP by choosing dt in the spectral representation. This is the
  ## latest version of the implementation that unties variables for the
  ## spectrum computation from those for generating the GP.

  ## cov - function handle to covariance funtion
  ## rng_seed - seed for random number generator

  ## cov=@(t) sigma^2*exp(-t.^2/(2*gamma^2));
  ## options = optimoptions('fsolve','TolFun',1e-8);
  ## ^^ set param 'TolFun'=1e-8 for the 'fsolve' optimizer

  ## compute spectrum and spectral cut-off frequency
  Nt_cor <- 100
  dt <- l / Nt_cor
  eps.cov <- 1e-4
  Tmax <- pracma::fzero(function(t) cov.fn(t, l = l) - eps.cov,
                        5 * l, tol = 1e-8)

  Tmax <- ceiling(Tmax$x / dt) * dt
  Nt <- Tmax / dt
  ind_0 <- Nt + 1
  k <- -Nt:(Nt - 1)
  tcov <- k * dt
  cov.t <- cov.fn(tcov,l=l)

  freq <- k / (2 * Tmax);
  w <- 2 * pi * freq;
  dw <- 2 * pi / (2 * Tmax);
  S <- abs(pracma::fftshift(stats::fft(cov.t)) * Tmax / Nt) / (2 * pi);
  Int_S <- cumsimpson(w,S);

  ## Test for convergence of spectral integral
  rel_error_Int_S = (Int_S[(2*Nt - 10L):(2*Nt)] -
                       Int_S[(2*Nt - 11L):(2*Nt - 1L)]) /
    Int_S[(2*Nt - 10L):(2*Nt)]
  eps_Int_S <- 1e-3
  if (sum(abs(rel_error_Int_S) > eps_Int_S) > 0) {
    stop("w-range for spectrum too small. Increase Tmax. (eps_Int_S = ",
         eps_Int_S, ").")
  }

  eps_cutoff <- 1e-3

  ind_wc <- which(Int_S > (1 - eps_cutoff) * Int_S[2*Nt])[1]
  Nw <- ind_wc - ind_0
  wc <- Nw * dw

  ## generate stochastic process with dt close to intial guess di_GP
  Nw <- max(ceiling((steps+1)*dt.GP*wc/(2*pi)),2)
  # # print(Nw)
  # Nw <- Nw.input
  # Nw = 128
  dw2 <- wc/Nw
  T0 <- 2 * pi / dw2
  Mw <- floor(T0 / dt.GP)
  dt <- T0 / Mw

  if (dt > pi / wc)
    stop('Reduce dt to be consistent with cut-off frequency wc')

  ## recompute spectrum with new spetracl discretisation
  dts <- 2 * pi / (2 * wc)
  cov_t <- cov.fn((-Nw:(Nw - 1)) * dts, l=l)
  S2 <- abs((stats::fft(cov_t) * dts)) / (2 * pi)
  w2 <- (-Nw:(Nw - 1)) * dw2

  A <- vector("numeric", Mw)
  A[1:Nw] <- sqrt(2 * S2[1:Nw] * dw2)
  A[1] <- 0

  ## not sure why the seed is being set? so I skipped
  # if (!is.null(rng_seed)) set.seed(rng_seed)

  phi <- vector("numeric", Mw)
  phi[1:Nw] <- 2 * pi * stats::runif(Nw)
  B <- sqrt(2) * A * exp(complex(real = 0, imaginary = 1) * phi)

  GP <- Mw * Re(stats::fft(B, inverse = TRUE) / length(B))
  return(list("GP"=GP[1:(steps+1)],"dt"=dt))
}


#' Function for Dirac Delta function.
#'
#' @param x x
#' @param y y
#'
#' @return Dirac Delta function
#' @export
#'
#' @examples
delta<-function(x=0,y=0){
  if (x==y){
    return(1)
  } else {
    return(0)
  }
}

#' Calculate the covarience of two given points.
#'
#' @param x  x
#' @param y  y
#' @param sigma_f2  signal variance
#' @param sigma_v2 noise variance
#' @param l length scale
#'
#' @return covariance
#' @export
#'
#' @examples
covariance <- function(x,y,sigma_f2,sigma_v2,l){
  cov <- sigma_f2 * exp( (-(x-y)**2)/(2* l**2) )   +sigma_v2* delta(x,y)
  return(cov)
}

#' Function that returns the covariance matrix Sigma.
#'
#' @param tStop end time
#' @param l length scale
#' @param sigma_f2 signal variance
#' @param sigma_v2 noise variance
#' @param steps steps
#'
#' @return Covariance matrix
#' @export
#'
#' @examples
GetSigma <- function(tStop,l,sigma_f2,sigma_v2,steps){
  grid <- seq(0,tStop,tStop/steps)
  Sigma <- matrix(c(1:((length(grid))**2)),nrow=(length(grid)))
  for (i in 1:(length(grid))){
    for(j in 1:(length(grid))){
      Sigma[i,j]<- covariance(grid[i],grid[j],sigma_f2,sigma_v2,l)
    }
  }
  return(Sigma)
}

# A function to calculate the pi(l|x,y) = pi(l|x) (l doesn't depend of y).
#' Title
#'
#' @param l Length scale
#' @param x Intensity function
#' @param end.time Experiment end time
#' @param sigma.f2 signal variance
#' @param sigma.n2 noise variance
#'
#' @return Marginal of the length scale
#' @export
#'
#' @examples
log_pi_l <- function(l,x,end.time,sigma.f2,sigma.n2){
  # Need to consider the logarithm of the intensity as prior is defined on the logarithm.
  x <- log(x)

  # First calculate the prior for l (l~exp(0.01))
  log.prior.l <- stats::dexp(l, rate = 0.01, log = TRUE)
  # log.prior.l <- dgamma(l, shape = 2.5, rate = 0.5, log = TRUE)
  # log.prior.l <- dlnorm(l, meanlog = log(8), sdlog = 0.3, log = TRUE)


  # Next need to calculate the likelihood of the intensity given l. pi(x|l).
  # First get the covariace matrix Sigma.l.
  Sigma.l <- GetSigma(end.time,l,sigma.f2,sigma.n2,(length(x)-1))
  # Calculate the density at x of the MVN(0,Sigma.l) distribution.
  log.likeli.l <-mvtnorm::dmvnorm(x, mean = rep(0,length(x)), sigma = Sigma.l, log = TRUE)

  return(log.prior.l + log.likeli.l)
}


#' Title za
#'
#' @param d.spikes za
#' @param end.time za
#' @param iter za
#' @param burn za
#' @param steps za
#' @param w za
#' @param hyper.param za
#' @param prior.x za
#' @param start.l za
#' @param ISI.type za
#' @param start.hyper za
#' @param hyper.start za
#' @param x.start za
#' @param sigma.h  za
#' @param sigma.l  za
#' @param T.min.param za
#' @param sigma.t  za
#' @param l.param za
#' @param initial.l  za
#'
#' @return za
#' @export
#'
#' @examples
mcmc_gp <- function(d.spikes, end.time, iter, burn, steps, w, prior.x =c(1000,0.001), hyper.param = c(1,0.001), ISI.type = "Gamma", sigma.h =NA, start.hyper = 1000, hyper.start = 1, l.param = 1, start.l = 1000, sigma.l = NA, initial.l = 1, T.min.param = NA, sigma.t = NA, x.start = NULL){
  # Error messages
  # Check hyper.param has length 1 or 2. (for either fixed at X or prior is Gam(X,X))
  if(length(hyper.param)>2){stop("Error! The length of hyper.param != 1 or 2.\n The hyper.param you entered has length ",
                                 length(hyper.param),". \n")}

  # Set up the prior parameters.
  sigma.f2 <- prior.x[1]
  sigma.n2 <- prior.x[2]

  out.likeli <- NULL
  out.piL <- NULL

  # Update to store the number of times that we accept the candidate value.
  update <- c(0,0)

  # Setup if we need to sample T.min.
  if(is.na(T.min.param)){
    T.min = 0 ;  max.T.min = 0.1
  }
  else if(length(T.min.param)==1 && is.numeric(T.min.param)){
     # Fixed T.min value
     T.min = T.min.param

     # Define T.max so that we don't put forward a value larger then possible.
     max.T.min <- 100000
     for(i in 1:ncol(d.spikes)){
       s <- d.spikes[,i]
       s <- s[!is.na(s)]
       new.T.max <- min(s[-1] - s[-length(s)])
       max.T.min <- min(max.T.min,new.T.max)
     }

     if(T.min > max.T.min){
       stop('Your T.min value is larger then the smallest ISI.')
     }
   }
  else if(length(T.min.param)==2 && is.numeric(T.min.param)){
    out.T.min <- rep(NA,iter)
    if(is.na(sigma.t) == TRUE){stop("Error! If you want to use mcmc to get T.min estimation you require a sigma.t value.")}

    # Define T.max so that we don't put forward a value larger then possible.
    max.T.min <- 100000
    for(i in 1:ncol(d.spikes)){
      s <- d.spikes[,i]
      s <- s[!is.na(s)]
      if(s[length(s)]-end.time <1e-10){s <- s[-length(s)]}  # Remove last time if it corresponds to length of experiment

      new.T.max <- min(s[-1] - s[-length(s)])
      max.T.min <- min(max.T.min,new.T.max)
    }

    T.min = min(1e-10, max.T.min/2) # Set T.min to a valid value.
  }
  else{
    stop('Invalid Refractory period settings!')
  }

  # Setup if we need to sample l.
  if(length(l.param) == 1 && is.numeric(l.param)){
    l.cur = l.param
  }
  else if(length(l.param) == 2 && is.numeric(l.param)){
    # Check valid sigma.l
    if(!(length(sigma.l) == 1 && is.numeric(sigma.l))){
      stop('Invalid sigma.l!')
    }
    if(sigma.l < 0){
      stop('sigma l must be positive!')
    }

    # Check valid start.l
    if(!(length(start.l) == 1 && is.numeric(start.l))){
      stop('Invalid start.l!')
    }

    l.cur = initial.l
    out.l <- rep(NA,iter)
  }

  ###############
  ## Set the initial value of x. (centered about the mean spiking intensity (#spikes/time), with a wave (sine))
  # wheere if x.start is defined this is the starting value (x.start must be a vector of length steps+1)
  ###############
  t_orig <- seq(0,end.time, end.time/steps)
  avg.spikes <- sum(!is.na(d.spikes))/(ncol(d.spikes)*end.time)
  # x.cur <- rep(avg.spikes,steps+1)
  if(is.null(x.start)){
    x.cur <- rep(avg.spikes,steps+1) + (0.5*avg.spikes)*sin((t_orig*5)/end.time)  #  Begin with wave so that the inital length scale is appropiate.
  }
  else{
    x.cur <- x.start
  }

  # Set the inital value of the hyper parameters if required.
  # I.e if hyper.param has length 1 set hyp.cur to this value. If not set hyp.cur to hyper.start value.
  if(length(hyper.param) ==1){
    hyp.cur <- hyper.param
  }
  else(hyp.cur <- hyper.start)


  # Create stores for the outputs, matrix out.x for x, and vectors out.l, out.hyper for l and ISI hyper parameter respectively.
  out.x <- matrix(NA, nrow=iter, ncol=steps+1)
  out.hyper <- rep(NA,iter)

  # When calculating l we need to shorten the length of the intensity function. We do this to 50.
  factor <- steps / 50
  shortened <- seq(1,(steps+1), factor)

  pb <- utils::txtProgressBar(min = 1, max = iter+burn, style = 3)
  # mcmc loop starts here
  for (i in 2:(iter+burn)) {

    ###############
    # Update x
    ###############
    # Draw from the 'prior' i.e \nu \sim \mathcal{N}(\mathbf{0},\mathbf{\Sigma})
    nu <- stat_gp(cov.fn, l = l.cur, dt.GP = end.time/steps, steps = steps)$GP

    # Get the candidate x.
    x.can.star <- (1-w**2)**0.5 * log(x.cur)  + w * nu
    x.can <- exp(x.can.star)
    # print(d.spikes)
    # print(T.min)
    # print(end.time)
    # print(x.can)
    # calculate log likelihood.
    log.likli.can <- log_pi_x(d.spikes=d.spikes, hyper.param = hyp.cur, x= x.can, end.time = end.time, T.min = T.min, max.T.min = max.T.min, ISI.type = ISI.type, do.log = TRUE)
    log.likli.cur <- log_pi_x(d.spikes=d.spikes, hyper.param = hyp.cur, x= x.cur, end.time = end.time, T.min = T.min, max.T.min = max.T.min, ISI.type = ISI.type, do.log = TRUE)
    # draw u from U(0,1).
    u <- stats::runif(1)

    # M-H ratio
    if (log(u) < log.likli.can - log.likli.cur) {
      update[1] <- update[1]+1
      x.cur <- x.can

    }





    ###############
    # update l (doesn't depend of ISI parameters)
    ###############
    # Check we input l as a parameter.
    if(length(l.param) == 2){
      # Only start updating l once we have reached the start.l-th iteration.
      # Only update every 2 (default) iterations change i%%2 == 0 if you want to update a different no of times.
      if(i > start.l && i%%2 == 0){

        # Get the candidate value for l.
        l.can <- stats::rnorm(1,l.cur,sigma.l)

        # Check the candidate value is in a reasonable range (must be >0 and < 50 chosen ad hoc, may need changing?)
        if (l.can >0 && l.can < 50){

          # Calculate the likelihood.
          # Note this is calculated with a shortened intensity to aid computation time.
          log.pi.l.can <-  log_pi_l(l.can,x.cur[shortened],end.time,sigma.f2,sigma.n2)
          log.pi.l.cur <-  log_pi_l(l.cur,x.cur[shortened],end.time,sigma.f2,sigma.n2)
          # draw u from U(0,1)
          u <- stats::runif(1)

          # M-H ratio
          if (log(u) < log.pi.l.can - log.pi.l.cur) {
            update[2] <- update[2]+1
            l.cur <- l.can
          }
        }
      }

    }




    ###############
    # Update hyper parameter.
    ###############
    if(length(hyper.param) == 2 && i > start.hyper){
      # Get the candidate value for the hyper.
      hyp.can <- stats::rnorm(1, hyp.cur, sigma.h)

      # Chec candidate value is reasonable.
      if (hyp.can > 0) {
        log.pi.cur <- log_pi_hyper_param(d.spikes = d.spikes, hyper.param = hyp.cur, x=x.cur, end.time = end.time, prior.hyper.param = hyper.param, T.min = T.min, max.T.min = max.T.min, ISI.type = ISI.type)
        log.pi.can <- log_pi_hyper_param(d.spikes = d.spikes, hyper.param = hyp.can, x=x.cur, end.time = end.time, prior.hyper.param = hyper.param, T.min = T.min, max.T.min = max.T.min, ISI.type = ISI.type)

        # M-H ratio

        # draw u from U(0,1)
        u <- stats::runif(1)

        if (log(u) < log.pi.can - log.pi.cur) {
          hyp.cur <- hyp.can
        }
      }
    }

    ###############
    ### Update T.min
    ###############
    if(length(T.min.param)==2){
      T.min.can <- stats::rnorm(1, T.min, sigma.t)
      if (T.min.can > 0 && T.min.can < max.T.min) {

         log.pi.T.cur <- log_pi_Tmin(d.spikes = d.spikes, T.min = T.min, max.T.min = max.T.min, hyper = hyp.cur, x = x.cur, end.time = end.time, prior.T.min = T.min.param, ISI.type =ISI.type)
        log.pi.T.can <- log_pi_Tmin(d.spikes = d.spikes, T.min = T.min.can, max.T.min = max.T.min, hyper = hyp.cur, x = x.cur, end.time = end.time, prior.T.min = T.min.param, ISI.type =ISI.type)
        # M-H ratio

        # draw from a U(0,1)
        u <- stats::runif(1)

        if (log(u) < log.pi.T.can - log.pi.T.cur) {
          T.min <- T.min.can
        }
      }
    }



    ###############
    # Store results.
    ###############

    if (i>burn){
      j <- i-burn
      out.x[j,]<- x.cur
      if(length(l.param)==2){out.l[j] <- l.cur}
      if(length(hyper.param)==2){ out.hyper[j] <- hyp.cur}
      if(length(T.min.param)==2){ out.T.min[j] <- T.min}

      out.likeli <- c(out.likeli,log.likli.cur)
      if(length(l.param)==2){out.piL <- c(out.piL, log.pi.l.cur )}
    }

    utils::setTxtProgressBar(pb, i) # update text progress bar after each iter
  }
  # Calculate the probability a move gets accepted.
  # print(update)
  p.x.acc = update[1]/(iter+burn)
  p.l.acc = update[2]/(iter+burn)

  hyp <- list(hyper = hyper.param) ; l <- list(l = l.param) ; t <- list(T.min = T.min.param)
  intensity <- list(intensity = out.x)
  if (length(hyper.param)==2){ hyp <- list(hyper = out.hyper)}
  if (length(T.min.param)==2){ t <- list(T.min = out.T.min)}
  if(length(l.param)==2){ l <- list(l = out.l)}

  output <- do.call("c",list(intensity,hyp,t,l))

  return(output)
}

