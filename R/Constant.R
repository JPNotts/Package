## Constant inference.

#' MCMC algorithm when the intensity function is constant.
#'
#' @param end.time end.time
#' @param prior.hyper prior.hyper
#' @param prior.x prior.x
#' @param iter iter
#' @param burn burn
#' @param ISI.type ISI.type
#' @param show.prog show.prog
#' @param prior.tmin prior.tmin
#' @param sigma.hyper sigma.hyper
#' @param sigma.x sigma.x
#' @param sigma.t.min sigma.t.min
#' @param d.spikes d.spikes
#'
#' @return mcmc
#' @export
#'
#' @examples
#' x <- 100
mcmc_constant <- function(d.spikes, end.time, prior.hyper = c(1,0.01), prior.x = c(1,0.01), prior.tmin = c(1,0.01), sigma.hyper, sigma.x, sigma.t.min, iter, burn, ISI.type, show.prog=T){
  # Create outputs and create initial values which stores the iterations.
  if(length(prior.x) == 1){
    out.x <- prior.x
    x.cur = prior.x
  }
  else if(length(prior.x) == 2){
    out.x <- rep(NA, iter)
    x.cur <- 1
  }
  else{
    stop('Invalid prior.x')
  }

  if(length(prior.hyper) == 1){
    out.hyper <- prior.hyper
    hyper.cur = prior.hyper
  }
  else if(length(prior.hyper) == 2){
    out.hyper <- rep(NA, iter)
    hyper.cur = 1
  }
  else{
    stop('Invalid prior.hyper')
  }

  if(length(prior.tmin) == 1){
    out.tmin <- prior.tmin
    tmin.cur = prior.tmin
  }
  else if(length(prior.tmin) == 2){
    out.tmin <- rep(NA, iter)
    tmin.cur = 1e-5 # Set tmin.cur to very small

  }
  else{
    stop('Invalid prior.tmin')
  }

  # Define T.max so that we don't put forward a value larger then possible.
  max.T.min <- 100000
  for(i in 1:ncol(d.spikes)){
    s <- d.spikes[,i]
    s <- s[!is.na(s)]
    if(s[length(s)]-end.time <1e-10){s <- s[-length(s)]}  # Remove last time if it corresponds to length of experiment
    new.T.max <- min(s[-1] - s[-length(s)])
    new.T.max
    max.T.min <- min(max.T.min,new.T.max)
  }


  if(show.prog){pb <- utils::txtProgressBar(min = 1, max = iter+burn, style = 3)}
  # Begin the mcmc loop.
  for (i in 1:(iter+burn)){

    ###########################
    ### x UPDATE
    ###########################
    if(length(prior.x) == 2){
      x.can <- stats::rnorm(1,x.cur, sigma.x)
      if(x.can > 0){

        log.pi.cur <- log_pi_x(d.spikes, hyper.cur, x= x.cur, end.time = end.time, T.min = tmin.cur, max.T.min = max.T.min, ISI.type =ISI.type, do.log = TRUE) + ((prior.x[1]-1) * base::log(x.cur) - prior.x[2]*x.cur)
        log.pi.can <- log_pi_x(d.spikes, hyper.cur, x = x.can, end.time = end.time, T.min = tmin.cur, max.T.min = max.T.min, ISI.type =ISI.type, do.log = TRUE) + ((prior.x[1]-1) * base::log(x.can) - prior.x[2]*x.can)

        # M-H ratio
        u <- stats::runif(1)

        if (log(u) < log.pi.can - log.pi.cur) {
          x.cur <- x.can
        }
      }
    }

    ###########################
    ### hyper UPDATE
    ###########################
    if(length(prior.hyper) == 2){
      hyper.can <-stats::rnorm(1,hyper.cur, sigma.hyper)
      if(hyper.can > 0){
        log.pi.cur <- log_pi_hyper_param(d.spikes, hyper.cur, x = x.cur, end.time = end.time, prior.hyper.param = prior.hyper, T.min = tmin.cur, max.T.min = max.T.min, ISI.type = ISI.type)
        log.pi.can <- log_pi_hyper_param(d.spikes, hyper.can, x = x.cur, end.time = end.time, prior.hyper.param = prior.hyper, T.min = tmin.cur, max.T.min = max.T.min, ISI.type = ISI.type)
        # M-H ratio
        u <- stats::runif(1)

        if (log(u) < log.pi.can - log.pi.cur) {
          hyper.cur <- hyper.can
        }
      }
    }

    ###########################
    ### tmin UPDATE
    ###########################
    if(length(prior.tmin) == 2){
      tmin.can <- stats::rnorm(1,tmin.cur,sigma.t.min)
      if(tmin.can >0 && tmin.can < max.T.min){
        log.pi.tmin.cur <- log_pi_Tmin(d.spikes, tmin.cur, max.T.min = max.T.min, hyper = hyper.cur, x=x.cur, end.time=end.time, prior.T.min = prior.tmin, ISI.type = ISI.type)
        log.pi.tmin.can <- log_pi_Tmin(d.spikes, tmin.can, max.T.min = max.T.min, hyper = hyper.cur, x=x.cur, end.time=end.time, prior.T.min = prior.tmin, ISI.type = ISI.type)

        # M-H ratio
        u <- stats::runif(1)

        if (log(u) < log.pi.tmin.can - log.pi.tmin.cur) {
          tmin.cur <- tmin.can
        }
      }
    }


    ###########################
    ### UPDATE iteration
    ###########################

    if ( i > burn){
      j <- i-burn
      if(length(prior.hyper) == 2){out.hyper[j] <- hyper.cur}
      if(length(prior.x) == 2){out.x[j] <- x.cur}
      if(length(prior.tmin) == 2){out.tmin[j] <- tmin.cur}
    }
   if(show.prog){ utils::setTxtProgressBar(pb, i) }# update text progress bar after each iter
  }


  hyp <- list(hyper = out.hyper) ; t <- list(T.min = out.tmin) ; intensity <- list(intensity = out.x)

  output <- do.call("c",list(intensity,hyp,t))

  return(output)
}
