#' Convert time series to spikes
#'
#' @param data data
#' @param time.step time.step
#' @param spike.occur spike.occur
#' @param val.change val change
#' @param reset reset
#'
#' @return spike times
#' @export
#'
#'
raw_spike_seq <- function(data, time.step, spike.occur, val.change = NULL, reset = spike.occur ){
  flag.spike <- FALSE
  spikes <- NULL
  spike.occ <- spike.occur[1] ; res <- reset[1]
  post.conc.max <- 1e5

  # Compute spikes if the given step is given by a single value.
  if(length(time.step) == 1){
    val.change <- c(val.change, length(data)*time.step +1)
    change.value <- val.change[1]
    counter <- 2
    for (index in 1:length(data)){

      #  Break loop of the experiment has finished (i.e first hit of NA)
      if( is.na(data[index]) ) {break}

      # Change the spike occur and reset values if we hit critical times.
      if( index * time.step >= change.value){
        spike.occ <- spike.occur[counter]
        res <- reset[counter]
        change.value <- val.change[counter]
        counter <- counter +1
      }

      # Reset the flag spike to false once the Ca concentration returns to below reset.
      if( flag.spike == TRUE && data[index] < res ){
        flag.spike <- FALSE
      }

      # Record the spike once the concentration hits spike.occur.
      if( flag.spike == FALSE && data[index] > spike.occ){
        spikes <- c(spikes,index * time.step)
        flag.spike <- TRUE
      }
    }
    # Adjoins the end time of the experiment.
    spikes <- c(spikes, length(data[!is.na(data)])*time.step)
  }

  # Compute spikes where we have a vector of length data of the time at each step.
  if(length(time.step) == length(data)){
    val.change <- c(val.change, time.step[length(time.step)] +1)
    change.value <- val.change[1]
    counter <- 2
    for (index in 1:length(data)){

      #  Break loop of the experiment has finished (i.e first hit of NA)
      if( is.na(data[index]) ) {break}

      # Change the spike occur and reset values if we hit critical times.
      if(time.step[index] >= change.value){
        spike.occ <- spike.occur[counter]
        res <- reset[counter]
        change.value <- val.change[counter]
        counter <- counter +1
      }

      # Reset the flag spike to false once the Ca concentration returns to below reset.
      # And add the spike time to the spikes (time at max thres since flag.spike=T)
      if( flag.spike == TRUE && data[index] < res ){
        spikes <- c(spikes,time.step[pos.conc.max])
        flag.spike <- FALSE
      }

      # Once threshold is hit save this index and set flag = T.
      if( flag.spike == FALSE && data[index] > spike.occ){
        # spikes <- c(spikes,time.step[index])
        pos.conc.max <- index
        flag.spike <- TRUE
      }

      if(flag.spike == TRUE && data[index] > data[pos.conc.max]){
        pos.conc.max <- index
      }

    }
    # Adjoins the end time of the experiment.
    spikes <- c(spikes, time.step[length(data[!is.na(data)])])
  }


  return(spikes)
}


#' Convert spike times to ISI times
#'
#' @param spikes spike times
#' @param adjust Flag, where if true we do not include the time to first spike
#'
#' @return ISI times
#' @export
#'
#' @examples
#' spikes <- c(3,9,12,15,17,21)
#' spike_to_ISI(spikes) #Keep time to first spike
#' spike_to_ISI(spikes,adjust = TRUE) #Remove time to first spike
spike_to_ISI <- function(spikes, adjust=FALSE){
  pre.spikes <- c(0,spikes[-length(spikes)])
  ISI <- spikes - pre.spikes

  if(adjust == TRUE){   # Adjust means remove first ISI (i.e start at time of first spike)
    ISI <- ISI[-1]
  }

  return(ISI)
}


#' Convert ISI times to spike times
#'
#' @param ISI  ISI times.
#'
#' @return Spike times.
#' @export
#'
#' @examples
#' ISI <- c(1,3,2,4,1,5)
#' ISI_to_spike(ISI)
ISI_to_spike <- function(ISI){
  return(cumsum(ISI))
}


#' Remove transient periods
#'
#' @param spikes spikes with end.time adjoined to the end
#' @param time the region/s of transients
#'
#' @return spikes with transients removed
#' @export
#'
#' @examples
#' spikes <- c(0.1,0.3,0.8,0.9,2,4,8,15,18,23)
#'
#' remove_transients(spikes, 1) # remove transients at start
#' remove_transients(spikes,c(0,7.5,16)) # remove transients in middle
#' remove_transients(spikes, c(0,16,25)) # remove transients at end
remove_transients <- function(spikes,time){
  # if(time[(length(time))] > spikes[length(spikes)]){
  #   stop("You're trying to cut off time after the end time. (end time = ",spikes[length(spikes)],")")}
  cuts.left <- (length(time)+1)/2

  # Cut the end section off if the interval contains the end time.
  if(time[length(time)] > spikes[length(spikes)] && length(time) > 2){
    new.end <- time[length(time) - 1]
    spikes <- spikes[which(spikes < new.end)]
    cuts.left <- cuts.left - 1

    # Add the new end.time?
    spikes <- c(spikes, new.end)
  }

  # Loop to cut all the sections not at the beginning.
  while(cuts.left >1.5){
    t.lower <- time[2*(cuts.left-1)] ; t.upper <- time[2*cuts.left-1]
    last.spike <- max(which(spikes < t.lower))
    next.spike <- min(which(spikes > t.upper))
    time.gap <- spikes[next.spike] - spikes[last.spike]
    spikes <- c(spikes[1:last.spike], spikes[(next.spike+1):length(spikes)]-time.gap )
    cuts.left <- cuts.left-1
  }

  # Remove the transients at the start, corresponds to time before time[1].
  if( time[1] > spikes[1]){        # Check the time occurs after the first spike
    last.spike <- max(which(spikes < time[1]))
    start.time <- spikes[last.spike+1]
    spikes <- spikes[-(1:(last.spike+1))]-start.time
  }

  return(spikes)
}


#' Remove Linear trends
#'
#' @param ISI The ISI times
#' @param rm.end.time Flag for whether we need to remove the last entry (end.time) of the ISI times.
#'
#' @return ISI.times with linear trend removed.
#' @export
#'
#' @examples
#' ISI <- c(1,3,2,4,1,5)
#' ISI.no <- 1:6
#' remove_trend(ISI, rm.end.time=FALSE)
#' stats::lm(ISI.no ~ remove_trend(ISI, rm.end.time=FALSE)) # Show the linear trend has been removed.
remove_trend <- function(ISI, rm.end.time = TRUE){
  spike.no <- 1:length(ISI)
  if(rm.end.time){
    if(length(ISI) < 3){return(ISI)}
    regression <- stats::lm(ISI[1:(length(ISI)-1)] ~ spike.no[1:(length(spike.no)-1)] ) # -1 to remove using the end time.
    linear.coeff <- regression$coefficients[2]
    normalise <- (spike.no-1)*linear.coeff
  }
  else{
    if(length(ISI) == 1){return(ISI)}
    regression <- stats::lm(ISI[1:(length(ISI))] ~ spike.no[1:(length(spike.no))] )
    linear.coeff <- regression$coefficients[2]
    normalise <- (spike.no)*linear.coeff
  }

  ISI.new <- ISI-normalise
  # Shift the ISI times to get the same mean.
  diff.sum <- sum(ISI.new)-sum(ISI)
  ISI.new <- ISI.new - diff.sum/length(ISI)
  # if(any(ISI.new<0)){
  #   print("Error! After removing trend you now have a negative ISI time!")
  # }
  return(ISI.new)
}


#' Merge spikes
#'
#' @param spikes  spike sequence
#' @param closeness Value which if two spikes are within this distance they are merged.
#' @param data time series data the spike sequence is derived from
#' @param time.step Vector of time from the time series data, or if step size is regular the time step of the experiment.
#'
#' @return Spike sequence where spikes within closeness have been merged.
#' @export
#'
#' @examples
spike_merge <- function(spikes, closeness = 40, data=NULL, time.step = NULL){
  if(length(spikes) <2){return(spikes)}
  # Case we merge spikes at their midpoint (No data provided)
  if(is.null(data)){
    for( i in (length(spikes)-1):1){
      if (spikes[i+1]-spikes[i]<closeness){
        spikes[i] <- (spikes[i]+spikes[i+1])/2
        spikes<- spikes[-(i+1)]
      }
    }
  }
  # Case we provide data and time.step. Merge by taking max of the conc of both spikes.
  else{
    if(length(time.step) == 1){
      time.step <- seq(time.step, length(data)*time.step, time.step)
    }

    for( i in (length(spikes)-2):1){ # Use  (i in (length(spikes)-1):1) to not exclude end time
      if (spikes[i+1]-spikes[i]<closeness){
        # Spikes too close.
        conc1 <- data[which(abs(time.step-spikes[i+1]) < 1e-5)]
        conc2 <- data[which(abs(time.step-spikes[i]) < 1e-5)]
        if(conc1 >= conc2){
          spikes <- spikes[-(i)]
        }
        else{
          spikes <- spikes[-(i+1)]
        }
      }
    }
  }

  return(spikes)
}


#' Get spikes from raw data
#'
#' @param data  time series data.
#' @param transients  transient periods where we want to remove spikes.
#' @param time.step Vector of time from the time series data, or if step size is regular the time step of the experiment.
#' @param spike.occur threshold where we record spikes.
#' @param val.change vector of times where the threshold value changes.
#' @param reset concentration at which we say a spike is over, default is set to spike.occur.
#' @param closeness Value which if two spikes are within this distance they are merged.
#' @param rm.trend Flag, where if TRUE we remove linear trends.
#'
#' @return spike sequence
#' @export
#'
#' @examples
data_to_spike <- function(data,transients = 2000, time.step = 2, spike.occur = 0.6,val.change =NULL, reset = spike.occur, closeness = 40, rm.trend = TRUE){

  spikes <- raw_spike_seq(data,time.step,spike.occur,val.change,reset)
  spikes <- remove_transients(spikes,transients)
  spikes <- spike_merge(spikes,closeness, data = data, time.step=time.step)
  if(rm.trend == TRUE){
    spikes <- ISI_to_spike( remove_trend( spike_to_ISI(spikes, adjust=FALSE) ) )
  }


  return(spikes)
}


#' Get spikes from .csv file of raw data.
#'
#' @param filename path to the csv file .
#' @param time.step Vector of time from the time series data, or if step size is regular the time step of the experiment.
#' @param transients  transient periods where we want to remove spikes.
#' @param spike.occur threshold where we record spikes.
#' @param closeness Value which if two spikes are within this distance they are merged.
#' @param rm.trend Flag, where if TRUE we remove linear trends.
#'
#' @return data frame containing spike sequences as columns.
#' @export
#'
#' @examples
csv_to_spike <- function(filename, time.step = 2, transients = 2000, spike.occur = 0.6, closeness = 40, rm.trend = TRUE){
  df <- utils::read.csv(filename)
  no.of.cell.data <- ncol(df)-1

  # Check all the variables are of the right form.
  if(length(transients) == 1){
    transients <- as.list(rep(transients,no.of.cell.data))
  }else if (length(transients) != no.of.cell.data){
    stop("You entered ",length(transients), " for the length of your transients vector and this does
       not equal the number of cell data ", no.of.cell.data, " or 1.")
  }
  if(length(spike.occur) == 1){
    spike.occur <- rep(spike.occur,no.of.cell.data)
  }else if (length(spike.occur) != no.of.cell.data){
    stop("You entered ",length(spike.occur), " for the length of your spike.occur vector and this does
       not equal the number of cell data ", no.of.cell.data, " or 1.")
  }
  # if(length(reset) == 1){
  #   reset <- rep(reset,no.of.cell.data)
  # }else if (length(reset) != no.of.cell.data){
  #   stop("You entered ",length(reset), " for the length of your spike.occur vector and this does
  #        not equal the number of cell data ", no.of.cell.data, " or 1.")
  # }


  for(i in 1:no.of.cell.data){
    print(i)
    # Split the spike.occur into spike occur and val.change.
    spike.occ <- spike.occur[[i]]
    leng <- length(spike.occ)
    if (leng > 1){
      val.change <- spike.occ[(leng/2+1.5):leng]
      spike.occ <-spike.occ[-((leng/2+1.5):leng)]
    }
    else{
      val.change <- NULL
    }


    # Convert data into spike for (i+1)th column.
    spikes <- data_to_spike(data = df[,i+1],
                            transients = transients[[i]],
                            time.step = time.step,
                            spike.occur = spike.occ,
                            val.change = val.change,
                            closeness = closeness,
                            rm.trend = rm.trend)
    # If on first iteration create data frame to store the spike times.
    if (i == 1){
      output <- data.frame(spikes)
    }else{

      # If the current spike seq is longer then any previous makke data frame same length by concaternating NA into data frame.
      if(length(spikes) > nrow(output)){
        output[ (nrow(output)+1):length(spikes),] <- NA
      }
      # If the current spike seq is shorter then any previous add NA to end of seq to make the same length.
      if(length(spikes) < nrow(output)){
        spikes <- c(spikes, rep(NA, (nrow(output)-length(spikes)) ))
      }

      # Add new column onto data frame for next spike sequence.
      output[ ,ncol(output)+1] <- spikes

    }


  }

  return(output)
}

