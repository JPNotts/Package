# QQ + KS plot stuff

#' A function to calclate q(t), the conditional intensity function.
#'
#' @param int.fn int.fn
#' @param end.time end.time
#' @param spike.times spike.times
#' @param ISI.type ISI.type
#' @param hyper hyper.param
#' @param t.min t.min
#' @param do.log do.log
#'
#' @return q(t)
#' @export
#'
#' @examples
#' 100
q.t <- function(int.fn, end.time, spike.times, ISI.type, hyper,t.min = 0, do.log = TRUE){
  steps <- length(int.fn)-1
  time.step <- end.time / steps

  # Calculate X, the integrated intensity function
  x <- inc.fineness(int.fn,4)
  X <- cumsum(x)*(end.time/(length(x)-1))

  # A vector to store the q.t.
  q <- rep(NA, steps + 1)

  spike.on.grid <- (spike.times%/%time.step)*time.step
  if(do.log == TRUE){
    for( i in 1:(steps+1) ) {

      t <- (i-1)*time.step
      if (any( abs(c(0,(spike.on.grid))-t)  < 1e-7 )){
        q[i] <- 0
        integral <- 0
        last.spike <- t
      }
      else{
        prob.cur <- exp(PDF_tmin(t,last.spike, hyper, end.time, x, X, ISI.type, t.min, do.log = T))
        # if(is.na(prob.cur)){stop('PDF nan')}
        integral <- integral + prob.cur * time.step
        q[i] <- prob.cur / ( 1 - integral)
         # cat('i = ',i,', prob.cur = ',prob.cur, ' and integral = ', integral,' giving q[i] =',q[i],'. \n')

      }

    }
  }
  else{
    for( i in 1:(steps+1) ) {
      # print(t)
      t <- (i-1)*time.step
      if (any( abs(c(0,(spike.on.grid))-t)  < 1e-7 )){
        q[i] <- 0
        integral <- 0
        last.spike <- t
      }
      else{
        prob.cur <- PDF_tmin(t,last.spike, hyper, end.time, x, X, ISI.type, t.min, do.log = F)
        integral <- integral + prob.cur * time.step
        # print(integral)
        q[i] <- prob.cur / ( 1 - integral)
      }

    }
  }

  return(q)
}

#' A function that returns R(t), the integrated intensity function.
#'
#' @param q q
#' @param time.step time.step
#'
#' @return
#' @export
#'
#' @examples
R.t <- function(q, time.step){
  return(cumsum(q)*time.step)
}

#' A function that returns the ISI times under the transformation.
#'
#' @param R.t R.t
#' @param spike.times spike.times
#' @param end.time end.time
#'
#' @return transformed ISI
#' @export
#'
#' @examples
get.tau <- function(R.t, spike.times,end.time){
  time.step <- end.time/(length(R.t)-1)
  spike.index <- spike.times / time.step

  tau <- rep(NA, length(spike.index)-1)
  for( i in 1:(length(spike.index)-1)){
    tau[i] <- R.t[spike.index[i+1]] - R.t[spike.index[i]]
  }

  return(tau)

}

#' A function that gets the qunatiles for the possion process of unit rate.
#'
#' @param N  N
#'
#' @return
#' @export
#'
#' @examples
model.quantiles <- function(N){
  s.k <- ( (1:N)-0.5 )/N
  t.k <- - log(1-s.k)
  return(t.k)
}


#' A function to increase the fineness of the intensity function, by a factor of the user's choice.
#'
#' @param int.fn int.fn
#' @param factor factor
#'
#' @return
#' @export
#'
#' @examples
inc.fineness <- function(int.fn, factor){
  new.length <- length(int.fn) + (length(int.fn)-1)*(factor-1)
  output <- rep(NA, new.length )

  new.pos <- seq(1,new.length,factor)
  output[new.pos]<- int.fn

  for ( i in 1:(length(int.fn)-1)){
    for (j in 1:(factor-1)){
      output[new.pos[i]+j] <- int.fn[i] + (j/factor) *(int.fn[i+1]-int.fn[i])
    }
  }

  return(output)
}


#' A function to get the quantiles given spikes and intensity function.
#'
#' @param int.fn int.fn
#' @param end.time end.time
#' @param spike.times spike.times
#' @param ISI.type ISI.type
#' @param hyper.param hyper.param
#' @param t.min t.min
#' @param do.log do.log
#'
#' @return experimental quantiles
#' @export
#'
#' @examples
exper.quantiles <- function(int.fn, end.time, spike.times, ISI.type, hyper.param, t.min = 0, do.log = FALSE){
  q <- q.t(int.fn = int.fn, end.time=end.time, spike.times=spike.times, ISI.type=ISI.type, hyper =hyper.param, t.min = t.min, do.log=do.log)
  q[which(q<0)] =0 #Note by definition q should be positive, therefor set numerical negative values to zero.
  time.step <- end.time/(length(int.fn)-1)
  R <- R.t(q, time.step)
  tau <- get.tau(R, spike.times, end.time)
  return(sort(tau))

}


#' Get the QQ and KS quantiles
#'
#' @param spikes.df spikes
#' @param int.fn int.fn
#' @param end.time end.time
#' @param ISI.type ISI.type
#' @param hyper.param hyper.param
#' @param t.min t.min
#' @param rm.last.spike rm.last.spike
#' @param do.log do.log
#'
#' @return quantiles
#' @export
#'
#' @examples
QQ_KS_data <- function(spikes.df, int.fn, end.time, ISI.type, hyper.param, t.min = 0, rm.last.spike = TRUE, do.log = TRUE){
  no.seq <- ncol(spikes.df)
  # Create vector to store the experimental quantiles, model quantiles, and s.k.
  all.s.k <- numeric()
  all.exper.quant <- numeric()
  all.model.quant <- numeric()

  all.s.k <- matrix(NA, ncol=no.seq, nrow=nrow(spikes.df)-1)
  all.exper.quant <- all.s.k
  all.model.quant <-  all.s.k


  # Iterate through all the spike sequences.
  for(i in 1:no.seq){
    # Get the current spike.
    spikes <- spikes.df[,i]
    spikes <- spikes[!is.na(spikes)]
    if(rm.last.spike == TRUE){
      spikes <- spikes[-length(spikes)]
    }
    N <- length(spikes) - 1

    exper.quant <- exper.quantiles(int.fn[i,], end.time[i],spikes, ISI.type, hyper.param[i],t.min, do.log = do.log)
    if(any(exper.quant < -1000)){browser()}
    all.exper.quant[1:length(exper.quant),i] <- exper.quant

    all.model.quant[1:N,i] <- model.quantiles(N)

    s.k <- ( (1:N)-0.5 )/N
    all.s.k[1:N,i] <-  s.k

    if(length(exper.quant) != length(s.k)) warning(paste0('Warning! length exper does not match model for seq ',i,'. Exper = ',length(exper.quant), 'and Model = ',length(s.k)))
    # print(i)
    # cat('exper length =', length(exper.quant), ' and model = ', length( model.quantiles(N)))
    # print('===============================')

  }

  output <- list(s.k = all.s.k, model = all.model.quant, exper = all.exper.quant)

  return(output)
  }

#' A function that returns the slopes of the individual spike cell data of QQ and KS plots.
#'
#' @param data data
#' @param spikes spikes
#'
#' @return slopes
#' @export
#'
#' @examples
get_slopes <- function(data, spikes){
  all.model.quant <- data$model
  all.exper.quant <- data$exper
  all.s.k <- data$s.k

  # Firstly we need to calculate pos to work out what values belongs to each cell.
  slopes.QQ <- NULL
  slopes.KS <- NULL
  for (i in 1:ncol(spikes)){
    current <- spikes[,i]

    # Get the correspoinding values from quantiles.
    # exper.cur <- all.exper.quant[start:(start+length.for.cell-1)]
    # model.cur <- all.model.quant[start:(start+length.for.cell-1)]
    # s.k.cur <- all.s.k[start:(start+length.for.cell-1)]

    exper.cur <- all.exper.quant[,i][!is.na(all.exper.quant[,i])]
    model.cur <- all.model.quant[,i][!is.na(all.model.quant[,i])]
    s.k.cur <- all.s.k[,i][!is.na(all.s.k[,i])]

    # Get QQ values.
    # Use lm to get the slope.
    regressionQQ <- stats::lm(model.cur ~ exper.cur)
    slopes.QQ <- c(slopes.QQ, regressionQQ$coefficients[[2]])

    # Get KS values.
    u.k.cur <- 1 - exp(-exper.cur)
    regressionKS <- stats::lm(s.k.cur ~ u.k.cur)
    slopes.KS <- c(slopes.KS, regressionKS$coefficients[[2]])

  }

  return(list(QQ = slopes.QQ, KS = slopes.KS))
}
