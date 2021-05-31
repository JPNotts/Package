# Functions used for analysis / plots
# Where we assume the data is assimilated into a single .RData file.


#' Maximum of a column
#'
#' @param data data
#'
#' @return Column maximum
#' @export
#'
#' @examples
colMax <- function(data){
  sapply(data, max, na.rm = TRUE)
}

#' Create rastor plot of the spikes
#'
#' @param spikes spikes
#' @param inc.end inc.end
#' @param add add
#' @param main main
#'
#' @return plot
#' @export
#'
rastor <- function(spikes, add = FALSE, main='', rm.empty = T){
  # Remove columns that have no spikes.
  has.spikes <- which(!is.na(spikes[1,]))
  if(rm.empty){
    spikes <- spikes[has.spikes]
  }

  # First create the plot
  if(!add){
    plot(0,0,type = "n", ylim = c(1,(ncol(spikes)+1)), xlim = c(0,max(spikes[!is.na(spikes)])),
         xlab = 'Time (s)',
         # xlab ='', xaxt='n',
         yaxt = "n", ylab = "", cex.lab = 1.8, cex.axis = 1.6, main=main)
  }

  # Add the rastor
  for(i in 1:ncol(spikes)){
    spikes.cur <- spikes[!is.na(spikes[,i]),i]
    spikes.cur <- spikes.cur[-length(spikes.cur)]  # Remove the end time.
    y.axis <- rep(i,length(spikes.cur))
    points(spikes.cur, y.axis, pch = 20,cex=1.3)
  }
  # add end of experiment time?
  end.times <- colMax(spikes)
  for(i in has.spikes){
    lines(c(end.times[i],end.times[i]), c(i-0.5,i+0.5), col=2 ,lwd=2)
  }

}


#' Check which files we have results for
#'
#' @param filepath
#'
#' @return
#' @export
#'
#' @examples
check_result <- function(filepath){
#Create a table to store which files we have.
spikes <- read.csv('spikes.csv')[,-1]
seq <- ncol(spikes)
rm(spikes)
ISI.types <- c('Gam','IG','LN','Wei','Poisson', 'Tmin')
models <- c('PWC', 'GP', 'Constant')
out <- matrix(NA, ncol=seq, nrow=length(models)*length(ISI.types))

for(j in 1:length(models)){
  files <- list.files(paste0(filepath,'/',models[j]))

  for(i in 1:length(ISI.types)){
    cur <-  files[grepl(ISI.types[i],files)]
    have.file <- sort(parse_number(cur))
    out[(j-1)*length(ISI.types) +i,have.file] = 1
  }
}
colnames(out) <- paste0('seq ',1:seq)
rownames(out) <- paste0(rep(models,each =length(ISI.types)),' ', rep(ISI.types,length(models)))

error.seq <-which(!apply(out,2,function(x) length(x[!is.na(x)]) %in% c(0,length(models)*length(ISI.types)) ))

return(list(all =out, issues_seq = error.seq))

}

#' Return QQ and KS for single ISI.type and model type
#'
#' @param filepath filepath
#' @param model model
#' @param ISI.type ISI.type
#'
#' @return QQ and KS
#' @export
#'
#' @examples
get_QQ_KS_singleISI <- function(filepath,model,ISI.type, min.spikes=0, ignore = NA,rescale=T, opthyp = T){
  # load in the results
  load(filepath)

  # Check ISI.type is valid
  if(!(ISI.type %in% c('Gamma', 'LogNormal', 'InverseGaussian', 'Exponential', 'Weibull'))){
    stop('Invalid ISI.type.')
  }
  # Check model type is valid, and we have the model in the data.
  if(!(model %in% c('PWC', 'GP', 'Constant'))){
    stop('Invalid model type')
  }
  if(!exists(toString(model))){ # check we have variable called 'PWC', 'GP' or 'Constant
    stop('Model type not provided in results!')
  }

  # Choose the data we want from the results
  data <- get(model)[[ISI.type]]

  # Do we ignore any of the spike sequences?
  vect <- rep(T,ncol(spikes))
  if(!is.na(ignore[1])) vect[ignore] <- F


  # Get the intensity functions (Need to remove NA rows, and only do this for spike sequences with more than min.spikes spikes.)
  if(model %in% c('GP', 'PWC')){
    required.cells <- stats::complete.cases(data$x_mean) & apply(spikes,2,function(x) length(x[!is.na(x)]) > min.spikes) &vect
    int.fns <- data$x_mean[required.cells,]
  }
  else{
    required.cells <- stats::complete.cases(data$x) & apply(spikes,2,function(x) length(x[!is.na(x)]) > min.spikes) &vect
    int.fns <- data$x[required.cells,8]
    # turn constant intensity into a vector.
    store <- NULL
    for(i in 1:length(int.fns)){
      store <- rbind(store, rep(int.fns[i],6001))
    }
    int.fns <- store

  }

  spikes.seq <- spikes[,required.cells]
  end.time <- colMax(spikes.seq)

  if(rescale==T){
    scale <- 20/end.time
    int.fns <- int.fns/scale
    spikes.seq <- spikes.seq*scale
    end.time=rep(20,length(end.time))
  }

  if(ISI.type != 'Exponential'){
    hyp <- data$ISI_param[required.cells,8]
    if(opthyp){
      # function to optimise
      likelihyp <- function(hyper, parameters ){
        spikes <- parameters[[1]] ; int.fn <- parameters[[2]] ; end.time <- parameters[[3]] ; ISI.type <- parameters[[4]]

        return(-log_pi_hyper_param(d.spikes = data.frame(spikes), hyper.param=hyper, x=int.fn, end.time = end.time,
                                   prior.hyper.param = c(1,0.01), T.min = NULL, max.T.min = NA, ISI.type = ISI.type))
      }

      ISI.orig <- c('Gamma', 'LN', 'Weibull', 'NewIG')
      ISI.new <- c("Gamma", 'LogNormal', 'Weibull', 'InverseGaussian')
      ISI <- ISI.orig[which(ISI.type == ISI.new)]

      hyp <- rep(NA,ncol(spikes.seq))
      for(kk in 1:ncol(spikes.seq)){
        s <- spikes.seq[,kk][!is.na(spikes.seq[,kk])] ; int.fn <- int.fns[kk,] ; end.time.cur <- max(s,na.rm=T)
        parameters <- list(s,int.fn,end.time.cur,ISI)

        hyp[kk] <- optim(par=1, fn =likelihyp, method='Brent', parameters = parameters, lower = 0, upper = 500)$par
      }
      print(hyp)
    }
  }


  ans <- QQ_KS_data(spikes.seq, int.fns, end.time, ISI.type = ISI.type, hyper.param = hyp,t.min = 0, rm.last.spike = TRUE, do.log = TRUE)


  return(list(quantiles = ans, req_spikes = required.cells))
}

#' Get QQ and KS details
#;
#' @param filepath filepath
#' @param model model
#' @param return_val return_val
#'
#' @return details
#' @export
#'
#' @examples
get_QQ_KS <- function(filepath, model,return_val = T, min.spikes = 0,ignore = NA,rescale=F, opthyp = F){

  # Get the data for each ISI type.
  gamma.data <- get_QQ_KS_singleISI(filepath,model,"Gamma", min.spikes = min.spikes, ignore=ignore,rescale = rescale, opthyp = opthyp)
  poisson.data <- get_QQ_KS_singleISI(filepath,model,"Exponential", min.spikes = min.spikes, ignore=ignore,rescale = rescale, opthyp = opthyp)
  invG.data <- get_QQ_KS_singleISI(filepath,model,"InverseGaussian", min.spikes = min.spikes, ignore=ignore,rescale = rescale, opthyp = opthyp)
  weibull.data <- get_QQ_KS_singleISI(filepath,model,"Weibull", min.spikes = min.spikes, ignore=ignore,rescale = rescale, opthyp = opthyp)
  LN.data <- get_QQ_KS_singleISI(filepath,model,"LogNormal", min.spikes = min.spikes, ignore=ignore,rescale = rescale, opthyp = opthyp)

  # get the spikes
  load(filepath)

  # Box plots
  # From QQ, KS plot data, calculate the slope for each line.
  A <- get_slopes(gamma.data$quantiles,spikes[gamma.data$req_spikes])
  B <- get_slopes(poisson.data$quantiles,spikes[poisson.data$req_spikes])
  C <- get_slopes(invG.data$quantiles,spikes[invG.data$req_spikes])
  D <- get_slopes(weibull.data$quantiles,spikes[weibull.data$req_spikes])
  E <- get_slopes(LN.data$quantiles,spikes[LN.data$req_spikes])

  # Data frames that contain the slopes for each ISI distribution
  d.QQ <- data.frame(Gam = A$QQ, Poi = B$QQ, InvG = C$QQ, Wei = D$QQ, LN = E$QQ)
  d.KS <- data.frame(Gam = A$KS, Poi = B$KS, InvG = C$KS, Wei = D$KS, LN = E$KS)

  if(return_val){
    QQ_KS <- list(gamma = gamma.data, exponential = poisson.data, inverseGaussian = invG.data, weibull = weibull.data, logNormal = LN.data)
    slopes <- list(QQ = d.QQ, KS = d.KS)
    return(list(QQ_KS = QQ_KS, slopes = slopes))
  }
  # pdf("Figure6-QQ.pdf")
  # mar.default <- c(5,4,4,2) + 0.1
  # par(mar = mar.default + c(0, 4, 0, 0))
  boxplot(d.QQ, ylab = TeX('Slope'), cex.lab = 1.2, cex.axis = 1.5,
          ylim = c(0,3),
          # ylim = c(0,max(d.QQ)),
          boxwex = 0.5, boxcol = "blue", whisklty = "solid", medcol = "red", medlwd = 1, main ='QQ')
  abline(h=1, lty = "dashed")
  # dev.off()

  # pdf("Figure6-KS.pdf")
  # mar.default <- c(5,4,4,2) + 0.1
  # par(mar = mar.default + c(0, 4, 0, 0))
  boxplot(d.KS, ylim = c(0,10), ylab = TeX('Slope'), cex.lab = 1.5,cex.axis = 1.5,
          boxwex = 0.5, boxcol = "blue", whisklty = "solid", medcol = "red", medlwd = 1, main = "KS" )
  abline(h=1, lty = "dashed")
  # dev.off()
}


#' Plot the important plots
#'
#' @param filepath filepath
#' @param model model
#'
#' @return save plots
#' @export
#'
#' @examples
plot_everything <- function(filepath, model){
  original.dir <- getwd()
  # load in the data
  load(filepath)

  # Check model type is valid, and we have the model in the data.
  if(!(model %in% c('PWC', 'GP', 'Constant', 'All'))){
    stop('Invalid model type')
  }
  if(model %in% c('PWC', 'GP', 'Constant')){
    if(!exists(toString(model))){ # check we have variable called 'PWC', 'GP' or 'Constant
      stop('Model type not provided in results!')
    }
    results <- get(model)
  }

  # Base directory the results.RData file is.
  base.dir <- paste0(substr(filepath,1,(nchar(filepath)-13)))
  new.dir <- paste0(base.dir,'plots')
  dir.create(new.dir, showWarnings = FALSE)
  new.dir <- paste0(new.dir,'/',model)
  dir.create(new.dir, showWarnings = FALSE)
  path <- paste0(new.dir,'/')
  setwd(path)

  # Plot the QQ and KS plots
  if(model %in% c('Constant','GP','PWC')){
    {
      dir.create('QQ_KS',showWarnings = FALSE)
      setwd('QQ_KS')
      QQ_KS <- get_QQ_KS(filepath,model)

      # QQ slopes
      grDevices::pdf('QQ_slopes.pdf')
      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 4, 0, 0))
      graphics::boxplot(QQ_KS$slopes$QQ, ylab = 'Slope', cex.lab = 1.2, cex.axis = 1.5,
                        # ylim = c(0,3),
                        ylim = c(0,max(QQ_KS$slopes$QQ)),
                        boxwex = 0.5, boxcol = "blue", whisklty = "solid", medcol = "red", medlwd = 1, main ='QQ')
      graphics::abline(h=1, lty = "dashed")
      grDevices::dev.off()

      # KS slopes
      grDevices::pdf('KS_slopes.pdf')
      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 4, 0, 0))
      graphics::boxplot(QQ_KS$slopes$KS, ylab = 'Slope', cex.lab = 1.2, cex.axis = 1.5,
                        # ylim = c(0,3),
                        ylim = c(0,max(QQ_KS$slopes$KS)),
                        boxwex = 0.5, boxcol = "blue", whisklty = "solid", medcol = "red", medlwd = 1, main ='QQ')
      graphics::abline(h=1, lty = "dashed")
      grDevices::dev.off()

      # plot the individual QQ and KS plots
      for(i in 1:length(QQ_KS$QQ_KS)){
        # QQ
        grDevices::pdf(paste0('QQ_',names(QQ_KS$QQ_KS)[i],'.pdf'))
        plot(QQ_KS$QQ_KS[[i]]$quantiles$model, QQ_KS$QQ_KS[[i]]$quantiles$exper, pch = 20, xlab = 'model', ylab = 'exper',
             cex.lab = 1.2, cex.axis = 1.5)
        lines(c(0,100), c(0,100), col= 2, lwd = 2)
        grDevices::dev.off()

        # KS
        u.k <- 1 - exp(-QQ_KS$QQ_KS[[i]]$quantiles$exper)
        grDevices::pdf(paste0('KS_',names(QQ_KS$QQ_KS)[i],'.pdf'))
        plot(QQ_KS$QQ_KS[[i]]$quantiles$s.k, u.k, pch = 20, xlab = 'model', ylab = 'exper',
             cex.lab = 1.2, cex.axis = 1.5)
        lines(c(0,100), c(0,100), col= 2, lwd = 2)
        grDevices::dev.off()
      }

      setwd(path)
    }
  }



  ## Plot all the intensities on one graph, with raw / stimulus
  if(model == 'All'){
    {

      cex.lab = 1.5 ; cex.axis = 1.2
      have_dat <- which(stats::complete.cases(Constant$Gamma$x)) #Assume we have data for all the same stuff.
      for(i in have_dat){
        spikes.cur <- spikes[,i]
        raw.cur <- raw[,i]
        stim.cur <- Stimulus[,i]

        # Plot the 3 graphs together
        j=1
        pdf(paste0('Summary_cell_',i,'_',names(PWC)[j],'.pdf'),width=7,height=14)
        op <- par(mfrow = c(3,1),
                    oma = c(5,6,0,0) + 0.1,
                    mar = c(0,0,0.2,1) + 0)

        # Plot the stimulus
        plot(recordingWindow[1:length(raw.cur[!is.na(raw.cur)])],stim.cur[!is.na(raw.cur)],type='l', ylab ='Stimulus',xaxt='n',xlab='',
             cex.lab = cex.lab,cex.axis=cex.axis )
        # Plot the raw with threshold
        plot(recordingWindow[1:length(raw.cur[!is.na(raw.cur)])],raw.cur[!is.na(raw.cur)],type='l', ylab ='Mean(ROI)',xaxt='n',xlab='',
             cex.lab = cex.lab,cex.axis=cex.axis )
        rug(spikes.cur, col=2,ticksize = 0.05, lwd=2)

        # Plot the posteriors.

         # ymax <- max(GP[[j]]$x_high[i,], PWC[[j]]$x_high[i,],Constant[[j]]$x[i,6])
         ymax = 0.04
        ymin <- 0
        end.time <- max(spikes.cur, na.rm = T)
        plot(0,0, type='n', xlab = 'Time(s)', ylab ='', ylim = c(0,ymax), xlim = c(0,end.time),  cex.lab = cex.lab, cex.axis=cex.axis)
        s.GP <- seq(0,end.time, end.time/(length(GP[[j]]$x_low[i,])-1))
        s.PWC <- seq(0,end.time, end.time/(length(PWC[[j]]$x_low[i,])-1))
        polygon(c(s.GP,rev(s.GP)), c(GP[[j]]$x_low[i,],rev(GP[[j]]$x_high[i,])), col=rgb(0,0,0,alpha=100,max=255), border=NA)
        polygon(c(s.PWC,rev(s.PWC)), c(PWC[[j]]$x_low[i,],rev(PWC[[j]]$x_high[i,])), col=rgb(255,0,0,alpha=100,max=255), border=NA)
        polygon(c(0,end.time,end.time,0), c(Constant[[j]]$x[i,2], Constant[[j]]$x[i,6], Constant[[j]]$x[i,6], Constant[[j]]$x[i,2]), col=rgb(0,0,255,alpha=100,max=255), border=NA )
        lines(s.GP, GP[[j]]$x_mean[i,], lwd=2, col = 1)
        lines(s.PWC, PWC[[j]]$x_mean[i,], lwd=2, col = 2 )
        lines(c(0,end.time), c(Constant[[j]]$x[i,8], Constant[[j]]$x[i,8]), lwd=2, col = 4)
        rug(spikes.cur, ticksize=0.05, lwd=2)
        mtext('Stim',side=2,line=3,outer=T,at=0.825)
        mtext('Ca conc',side=2,line=3,outer=T,at=0.5)
        mtext('Result',side=2,line=3,outer=T,at=0.125)
        mtext('Time(s)', side=1,line=3,outer=T)
        dev.off()
      }
      setwd(path)

    }
  }




}
