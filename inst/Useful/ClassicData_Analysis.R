# ClassicData_Analysis
## Put all the required data into a single file 
 # directory <- '/Users/jakepowell/Desktop/Tom_Data_Raw/Dynamic Stimulation Data/Waves/141204_HEK_CAR_sineP1200_3/'
 directory <- getwd()
 source('~/Package/R/Analysis.R')
 check <- check_result(directory)
 # check$issues_seq
     # View(check$all)
 ignore <- as.numeric(check$issues_seq)
 enough.spikes <- which(apply(spikes, 2, function(x){length(x[!is.na(x)])}) > 18)
 spikes.for.QQ.KS <- enough.spikes[which(!(enough.spikes %in% ignore))]
 {

  # 1) Collect all files into one .RData file. 
  rm(list = ls())
  spikes <- read.csv('spikes.csv')[,-1]
  raw <- read.csv('data.csv')
  recordingWindow <- raw[,1]
  raw <- raw[,-1]
  threshold <- read.csv('store.csv')[,-1]
  if(file.exists('stim.csv')){
    Stimulus <- read.csv('stim.csv')[,-1]
  }
  
  # Store the pwc results in one variable, PWC. 
  nspikes <- ncol(spikes)
  ISI.list <- c('Gamma', 'InverseGaussian','LogNormal','Exponential','Weibull', 'Exponential_Tmin')
  PWC <- list('Gamma'=NULL, 'InverseGaussian'=NULL,'LogNormal'=NULL,'Exponential'=NULL,'Weibull'=NULL, 'Exponential_Tmin' = NULL)
  ISI <-  c('Gamma', 'NewIG','LN','Poisson','Weibull', 'PoiTmin')
  for(j in 1:length(ISI)){
    # Fill 4 tables: x_mean, x_low,x_high, ISIparam
    x_mean <- matrix(NA, ncol=2001, nrow = nspikes)
    x_low <- matrix(NA, ncol=2001, nrow = nspikes)
    x_high <- matrix(NA, ncol=2001, nrow = nspikes)
    ISI_param <- matrix(NA,ncol=8,nrow = nspikes)
    
    for(i in 1:nspikes){
      # Fill 4 tables: x_mean, x_low,x_high, ISIparam
      if(file.exists(paste0('PWC/Index',ISI[j],i,'.Rdata'))){
        load(paste0('PWC/Index',ISI[j],i,'.Rdata'))
        x_mean[i,] <- x.mean$mean
        x_low[i,] <- x.mean$lower
        x_high[i,] <- x.mean$upper
        ISI_param[i,] <- hyper.result 
        if(ISI[j] == 'PoiTmin'){
          ISI_param[i,] = tmin.result
        }
      }
    }
    cur <- list(x_mean = x_mean, x_low=x_low,x_high = x_high, ISI_param=ISI_param)
    PWC[[j]] <- cur 
  }
  
  # store the Constant results in one variable, Constant.
  nspikes <- ncol(spikes)
  ISI.list <- c('Gamma', 'InverseGaussian','LogNormal','Exponential','Weibull')
  Constant <- list('Gamma'=NULL, 'InverseGaussian'=NULL,'LogNormal'=NULL,'Exponential'=NULL,'Weibull'=NULL,'Exponential_Tmin' = NULL)
  ISI <-  c('Gamma', 'NewIG','LN','Poisson','Weibull','PoiTmin')
  for(j in 1:length(ISI)){
    # Fill 4 tables: x_mean, x_low,x_high, ISIparam
    x <- matrix(NA, ncol=8, nrow = nspikes)
    ISI_param <- matrix(NA,ncol=8,nrow = nspikes)
    
    for(i in 1:nspikes){
      # Fill 4 tables: x_mean, x_low,x_high, ISIparam
      if(file.exists(paste0('Constant/Cons',ISI[j],i,'.Rdata'))){
        load(paste0('Constant/Cons',ISI[j],i,'.Rdata'))
        x[i,] <- x.result
        ISI_param[i,] <- hyper.result 
        if(ISI[j] == 'PoiTmin'){
          ISI_param[i,] = tmin.result
        }
      }
    }
    cur <- list(x = x, ISI_param=ISI_param)
    Constant[[j]] <- cur 
  }
  
  # Store the GP results in one variable, GP. 
  GP <- NULL
  nspikes <- ncol(spikes)
  ISI.list <- c('Gamma', 'InverseGaussian','LogNormal','Exponential','Weibull')
  GP <- list('Gamma'=NULL, 'InverseGaussian'=NULL,'LogNormal'=NULL,'Exponential'=NULL,'Weibull'=NULL,'Exponential_Tmin' = NULL)
  ISI <-  c('Gamma', 'NewIG','LN','Poisson','Weibull','PoiTmin')
  for(j in 1:length(ISI)){
    # Fill 4 tables: x_mean, x_low,x_high, ISIparam
    x_mean <- matrix(NA, ncol=2001, nrow = nspikes)
    x_low <- matrix(NA, ncol=2001, nrow = nspikes)
    x_high <- matrix(NA, ncol=2001, nrow = nspikes)
    ISI_param <- matrix(NA,ncol=8,nrow = nspikes)

    for(i in 1:nspikes){
      # Fill 4 tables: x_mean, x_low,x_high, ISIparam
      if(file.exists(paste0('GP/Index',ISI[j],i,'.Rdata'))){
        load(paste0('GP/Index',ISI[j],i,'.Rdata'))
        x_mean[i,] <- x.mean$mean  ; x_mean[i,2001] <- x_mean[i,2000]
        x_low[i,] <- x.mean$lower ; x_low[i,2001] <- x_low[i,2000]
        x_high[i,] <- x.mean$upper ; x_high[i,2001] <- x_high[i,2000]
        ISI_param[i,] <- hyper.result
        if(ISI[j] == 'PoiTmin'){
          ISI_param[i,] = tmin.result
        }
      }
    }
    cur <- list(x_mean = x_mean, x_low=x_low,x_high = x_high, ISI_param=ISI_param)
    GP[[j]] <- cur
  }
  
  save(recordingWindow,spikes,raw,threshold,PWC,Constant, GP, file='Results.Rdata')
  # save(spikes,PWC,Constant, GP, file='Results.Rdata')
  
  if(file.exists('stim.csv')){
    save(recordingWindow,spikes,raw,threshold,PWC,Constant, GP,Stimulus, file='Results.Rdata')
  }

  
}

## Then use single file to: 
## Rastor plots.
 {
   directory <- getwd()
   filepath <- paste0(directory,'/Results.Rdata')
   load(filepath)
   source('~/Package/R/Analysis.R')
   check <- check_result(directory)
   ignore <- as.numeric(check$issues_seq)
   enough.spikes <- which(apply(spikes, 2, function(x){length(x[!is.na(x)])}) > 18)
   spikes.for.QQ.KS <- enough.spikes[which(!(enough.spikes %in% ignore))]

   
   ## - Rastor plot of all spikes. 
   source('~/Package/R/Analysis.R')
   pdf('Rastor_of_spikes_all.pdf')
   rastor(spikes, rm.empty = T)
   dev.off()
   
   pdf('Rastor_of_spikes_QQ_KS.pdf')
   rastor(spikes[,spikes.for.QQ.KS], rm.empty = T)
   dev.off()  
 }

## - Create .csv file with all the ISI param details. 
{
  # Get hyper for the spikes used in QQ/KS.
  ISI.types <- c("Gamma", "InverseGaussian","LogNormal", "Weibull")
  models <- c('Constant', 'PWC', 'GP')
  all <- NULL
  for(j in 1:length(models)){
    model <- get(models[j])
    for(i in 1:length(ISI.types)){
      cur <- cbind(as.numeric(spikes.for.QQ.KS), rep(models[j],length(spikes.for.QQ.KS)), rep(ISI.types[i],length(spikes.for.QQ.KS)), model[[ISI.types[i]]]$ISI_param[spikes.for.QQ.KS,])
      
      all <-rbind(all,cur)
    }
  }
  
  all[,4:11] <- as.numeric(all[,4:11])
  data <- all
  regions <- paste0('[',round(as.numeric(data[,5]),digits=2),', ', round(as.numeric(data[,9]),digits=2),']')
  new <- cbind(data[,1:3],round(as.numeric(data[,11]),digits=4), round(as.numeric(data[,7]),digits=4), regions)
  
  colnames(new) <- c('Seq.Number', 'Prior', 'ISI.type', 'mean', 'median', '95% interval')
  head(new)
  write.csv(new, file='ISIparam(QQ_KS)_database.csv')
  
  
  # Get standard database with all details
  all.hyp <- NULL 
  files <- list.files(getwd(), recursive = T)
  for(i in 1:length(files)){
    file.cur <- files[i]
    file.cur
    if(grepl('Rdata',file.cur, fixed=T)){
      load(file.cur) 
      
      if(exists('hyper.result')){
        # Extract ISI dist.
        types <- c('Gamma', 'LN',  'IG', 'Wei')
        for(j in 1:length(types)){
          if(grepl(types[j],file.cur, fixed=T)){
            ISI.cur <- types[j]
          }
        }
        ISI.cur
        # Extract sequence number
        seq.no <- readr::parse_number(file.cur)
        
        # Extract method
        types <- c('Con', 'GP', 'PWC')
        for(i in 1:length(types)){
          if(grepl(types[i],file.cur, fixed=T)){
            meth.cur <- types[i]
          }
        }
        meth.cur
        
        hyp.cur <- c(seq.no,meth.cur,ISI.cur, hyper.result)
        all.hyp <- rbind(all.hyp,hyp.cur)
        rm(hyper.result)
      }
      
    }
  }
  
  colnames(all.hyp) = c('Seq Number', 'Prior', 'ISI,type','min','2.5%', '25%','50%', '75%', '97.5%', 'max', 'Mean')
  all.hyp[,4:11] <- as.numeric(all.hyp[,4:11])
  # all.hyp[all.hyp[,1] == 3,]
  data <- all.hyp
  
  regions <- paste0('[',round(data[,5],digits=2),', ', round(data[,9],digits=2),']')
  spikes <- read.csv('spikes.csv')[,-1]
  no.spikes <- apply(spikes,2, function(x) length(x[!is.na(x)])-1)
  
  new <- cbind(data[,1:3],round(data[,11],digits=4), round(data[,7],digits=4), regions, no.spikes[data[,1]])

  colnames(new) <- c('Seq.Number', 'Prior', 'ISI.type', 'mean', 'median', '95% interval', 'No. Spikes')
  head(new)
  write.csv(new, file='ISIparam_database.csv')
  
}

## - Plot all individual intensity plots (table for constant).
{
  dir.create(paste0(directory, '/Individual_intensity_functions'), showWarnings = F)
  
  # For Constant
  ISI.type <- names(Constant)
  Constant.database <- NULL
  for(i in 1:length(Constant)){
    for(j in 1:ncol(spikes)){
      if(!is.na(Constant[[i]]$x[j,1])){
        cur <- c(round(Constant[[i]]$x[j,8],digits=4), paste0('[',round(Constant[[i]]$x[j,2],digits=4),', ',round(Constant[[i]]$x[j,7],digits=4),']'),
                 round(Constant[[i]]$ISI_param[j,8],digits=4), paste0('[',round(Constant[[i]]$ISI_param[j,2],digits=4),', ',round(Constant[[i]]$ISI_param[j,7],digits=4),']'))
        Constant.database <- rbind(Constant.database, c(j,ISI.type[i],cur))
      }
    }
  }
  colnames(Constant.database) <- c('Spike seq', 'ISI type', 'mean intensity', '95% intensity', 'mean ISI param', '95% ISI param')
  write.csv(Constant.database, file = paste0(directory,'/Constant_database.csv'),row.names=FALSE)
  
  # For PWC model
  new.dir <- paste0(directory, '/Individual_intensity_functions/PWC')
  dir.create(new.dir, showWarnings = F)
  ISI.type <- names(PWC)
  for(i in 1:length(PWC)){
    for(j in 1:ncol(spikes)){
      if(!is.na(PWC[[i]]$x_mean[j,1])){ ## Make sure we have an inferred intensity
        
        # Get spikes and end time.
        s <- spikes[,j][!is.na(spikes[,j])]
        end.time <- max(s) ; s <- s[-length(s)]
        print(i); print(j); print(s) ; print(end.time)
        pdf(paste0(new.dir,'/',ISI.type[i],j,'.pdf'))
        mar.default <- c(5,4,4,2) + 0.1
        par(mar = mar.default + c(0, 1, 0, 0))
        t <- seq(0,end.time,end.time/2000)
        plot(t, PWC[[i]]$x_high[j,], type='n', xlab= 'Time (s)', ylab = 'Intensity', cex.lab = 2 , cex.axis = 1.8, ylim = c(0,max(PWC[[i]]$x_high[j,])))
        polygon(c(t,rev(t)), c(PWC[[i]]$x_low[j,],rev(PWC[[i]]$x_high[j,])), col = rgb(1,0,0,alpha=0.5), border=NA )
        lines(t, PWC[[i]]$x_mean[j,],col=rgb(0,0,0,alpha=0.9),lwd=4)
        rug(s,ticksize=0.05,lwd=2)
        dev.off()
      }
    }
  }
  
  # For GP model
  new.dir <- paste0(directory, '/Individual_intensity_functions/GP')
  dir.create(new.dir, showWarnings = F)
  ISI.type <- names(GP)
  for(i in 1:length(GP)){
    for(j in 1:ncol(spikes)){
      if(!is.na(GP[[i]]$x_mean[j,1])){ ## Make sure we have an inferred intensity
        
        # Get spikes and end time.
        s <- spikes[,j][!is.na(spikes[,j])]
        end.time <- max(s) ; s <- s[-length(s)]
        
        pdf(paste0(new.dir,'/',ISI.type[i],j,'.pdf'))
        mar.default <- c(5,4,4,2) + 0.1
        par(mar = mar.default + c(0, 1, 0, 0))
        t <- seq(0,end.time,end.time/2000)
        plot(t, GP[[i]]$x_high[j,], type='n', xlab= 'Time (s)', ylab = 'Intensity', cex.lab = 2 , cex.axis = 1.8, ylim = c(0,max(GP[[i]]$x_high[j,])))
        polygon(c(t,rev(t)), c(GP[[i]]$x_low[j,],rev(GP[[i]]$x_high[j,])), col = rgb(1,0,0,alpha=0.5), border=NA )
        lines(t, GP[[i]]$x_mean[j,],col=rgb(0,0,0,alpha=0.9),lwd=4)
        rug(s,ticksize=0.05,lwd=2)
        dev.off()
      }
    }
  }
}

 ## - Plot all individual intensity plots split with raw data / stimulus.
 {
   dir.create(paste0(directory, '/Individual_intensity_functions_with_data'), showWarnings = F)
   new.dir <- paste0(directory, '/Individual_intensity_functions_with_data/PWC')
   dir.create(new.dir, showWarnings = F)
   new.dir <- paste0(directory, '/Individual_intensity_functions_with_data/GP')
   dir.create(new.dir, showWarnings = F)
   new.dir <- paste0(directory, '/Individual_intensity_functions_with_data/Constant')
   dir.create(new.dir, showWarnings = F)
   dir.cur <- paste0(directory, '/Individual_intensity_functions_with_data')
   
   model = c('Constant', 'PWC', 'GP')
   ISI.type <- c( "Gamma" , "InverseGaussian",  "LogNormal", "Exponential" ,  "Weibull", "Exponential_Tmin")
   # Loop over all spikes.
   for(i in 1:ncol(spikes)){
     for(ii in 1:length(model)){
       for(iii in 1:length(ISI.type)){
         cat('spikes =',i, ', model = ', model[ii], ', ISI = ', ISI.type[iii],'. \n')
         # Get the intensity values.
         if(model[ii] =='Constant'){
           int.mean <- rep(get(model[ii])[[iii]]$x[i,8],2001)
           high <- rep(get(model[ii])[[iii]]$x[i,6],2001)
           low <- rep(get(model[ii])[[iii]]$x[i,2],2001)
         }else{
           int.mean <- get(model[ii])[[iii]]$x_mean[i,]
           high <- get(model[ii])[[iii]]$x_high[i,]
           low <- get(model[ii])[[iii]]$x_low[i,]
         }
         # If we don't have a fit skip 
         if(is.na(int.mean[1] )){next}
         
         # Get the spikes
         s <- spikes[,i][!is.na(spikes[,i])]
         end.time <- max(s) ; s <- s[-length(s)]
         
         # Get the raw data / stim data
         conc <- raw[,i]
         store <- threshold[i,]
         if(exists('Stimulus')){
           stim <- Stimulus[,i]
           
           # plot the intensity with data. 
           pdf(paste0(dir.cur,'/',model[ii],'/',ISI.type[iii],'_',i,'.pdf'))
           op <- par(mfrow = c(3,1),
                     oma = c(5,6,0,0) + 0.1,
                     mar = c(0,0,0.2,1) + 0)
           cex.lab = 1.5 ; cex.axis=1.3
           # Plot the raw with threshold
           plot(recordingWindow[1:length(conc[!is.na(conc)])],stim[!is.na(conc)],type='l', ylab ='stim',xaxt='n',xlab='',
                cex.lab = cex.lab,cex.axis=cex.axis )
           plot(recordingWindow[1:length(conc[!is.na(conc)])],conc[!is.na(conc)],type='l', ylab ='Mean(ROI)',xaxt='n',xlab='',
                cex.lab = cex.lab,cex.axis=cex.axis )
           rug(s,ticksize = 0.05,col=2)
           # Plot the posterior intensity function.
           grid <-seq(0,end.time,end.time/(length(low)-1))
           plot(grid,int.mean, ylim = c(min(low),max(high)), type='n',
                ylab='Intensity (spikes/s)',xlab = 'Time (s)',cex.lab= cex.lab, cex.axis=cex.axis)
           polygon(c(grid,rev(grid)), c(low,rev(high)), col=rgb(0,0,0,alpha=100,max=255), border=NA)
           lines(grid,int.mean)
           rug(s,ticksize = 0.05,col=2)
           mtext('Stimulus',side=2,line=3,outer=T,at=0.85)
           mtext('Raw',side=2,line=3,outer=T,at=0.5)
           mtext('Intensity',side=2,line=3,outer=T,at=0.15)
           mtext('Time(s)', side=1,line=3,outer=T)
           dev.off()
         }
         else{
           # plot the intensity with data. 
           pdf(paste0(dir.cur,'/',model[ii],'/',ISI.type[iii],'_',i,'.pdf'))
           op <- par(mfrow = c(2,1),
                     oma = c(5,6,0,0) + 0.1,
                     mar = c(0,0,0.2,1) + 0)
           cex.lab = 1.5 ; cex.axis=1.3
           # Plot the raw with threshold
           plot(recordingWindow[1:length(conc[!is.na(conc)])],conc[!is.na(conc)],type='l', ylab ='Mean(ROI)',xaxt='n',xlab='',
                cex.lab = cex.lab,cex.axis=cex.axis )
           rug(s,ticksize = 0.05,col=2)
           # Plot the posterior intensity function.
           grid <-seq(0,end.time,end.time/(length(low)-1))
           plot(grid,int.mean, ylim = c(min(low),max(high)), type='n',
                ylab='Intensity (spikes/s)',xlab = 'Time (s)',cex.lab= cex.lab, cex.axis=cex.axis)
           polygon(c(grid,rev(grid)), c(low,rev(high)), col=rgb(0,0,0,alpha=100,max=255), border=NA)
           lines(grid,int.mean)
           rug(s,ticksize = 0.05,col=2)
           mtext('Raw',side=2,line=3,outer=T,at=0.75)
           mtext('Intensity',side=2,line=3,outer=T,at=0.25)
           mtext('Time(s)', side=1,line=3,outer=T)
           dev.off()
         }
         
         
         
       }
     }
   }
   
 }
 
## - Plot ISI intensity comparisons.
{
  dir.create(paste0(directory, '/ISI_intensity_comparisions'), showWarnings = F)
  new.dir <- paste0(directory, '/ISI_intensity_comparisions/PWC')
  dir.create(new.dir, showWarnings = F)
  ISI.type <- names(PWC)
  for(i in 1:ncol(spikes)){
    # If we have spikes do:
    if(!is.na(spikes[1,i])){
      # Work out the maximum height required in the plot. 
      max.intensity <- -10
      for(j in 1:length(ISI.type)){
        max.intensity <- max(max.intensity, PWC[[j]]$x_high[i,], na.rm = T)
      }
      
      # If we don't have any intensities break out of the current loop. 
      if(max.intensity<0 || is.na(max.intensity)){
        next
      }
      
      # Get spikes and end time.
      s <- spikes[,i][!is.na(spikes[,i])]
      end.time <- max(s) ; s <- s[-length(s)]
      
      colo_shade <-viridis::viridis(length(ISI.type),alpha=0.5)  
      colo <-viridis::viridis(length(ISI.type),alpha=1)  
      
      pdf(paste0(new.dir,'/spike_seq_',i,'.pdf'))
      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 1, 0, 0))
      t <- seq(0,end.time,end.time/2000)
      plot(0,0, type='n', xlab= 'Time (s)', ylab = 'Intensity', cex.lab = 2 , cex.axis = 1.8, ylim = c(0,max.intensity), xlim = c(0,end.time))
      for(j in 1:length(ISI.type)){
        polygon(c(t,rev(t)), c(PWC[[j]]$x_low[i,],rev(PWC[[j]]$x_high[i,])), col = colo_shade[j], border=NA )
      }
      for(j in 1:length(ISI.type)){
        lines(t, PWC[[j]]$x_mean[i,],col=colo[j],lwd=4)
      }
      rug(s,ticksize=0.05,lwd=2)
      legend('topright', col = colo, legend = ISI.type, lty = rep(1,length(ISI.type)))
      dev.off()
      
    }
    
  }
  
  new.dir <- paste0(directory, '/ISI_intensity_comparisions/GP')
  dir.create(new.dir, showWarnings = F)
  ISI.type <- names(GP)
  for(i in 1:ncol(spikes)){
    # If we have spikes do:
    if(!is.na(spikes[1,i])){
      # Work out the maximum height required in the plot. 
      max.intensity <- -10
      for(j in 1:length(ISI.type)){
        max.intensity <- max(max.intensity, GP[[j]]$x_high[i,], na.rm = T)
      }
      
      # If we don't have any intensities break out of the current loop. 
      if(max.intensity<0 || is.na(max.intensity)){
        next
      }
      
      # Get spikes and end time.
      s <- spikes[,i][!is.na(spikes[,i])]
      end.time <- max(s) ; s <- s[-length(s)]
      
      colo_shade <-viridis::viridis(length(ISI.type),alpha=0.5)  
      colo <-viridis::viridis(length(ISI.type),alpha=1)  
      
      pdf(paste0(new.dir,'/spike_seq_',i,'.pdf'))
      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 1, 0, 0))
      t <- seq(0,end.time,end.time/2000)
      plot(0,0, type='n', xlab= 'Time (s)', ylab = 'Intensity', cex.lab = 2 , cex.axis = 1.8, ylim = c(0,max.intensity), xlim = c(0,end.time))
      for(j in 1:length(ISI.type)){
        polygon(c(t,rev(t)), c(GP[[j]]$x_low[i,],rev(GP[[j]]$x_high[i,])), col = colo_shade[j], border=NA )
      }
      for(j in 1:length(ISI.type)){
        lines(t, GP[[j]]$x_mean[i,],col=colo[j],lwd=4)
      }
      rug(s,ticksize=0.05,lwd=2)
      legend('topright', col = colo, legend = ISI.type, lty = rep(1,length(ISI.type)))
      dev.off()
      
    }
    
  }
  
  new.dir <- paste0(directory, '/ISI_intensity_comparisions/Constant')
  dir.create(new.dir, showWarnings = F)
  ISI.type <- names(Constant)
  for(i in 1:ncol(spikes)){
    # If we have spikes do:
    if(!is.na(spikes[1,i])){
      # Work out the maximum height required in the plot. 
      max.intensity <- -10
      for(j in 1:length(ISI.type)){
        max.intensity <- max(max.intensity, GP[[j]]$x_high[i,], na.rm = T)
      }
      
      # If we don't have any intensities break out of the current loop. 
      if(max.intensity<0 || is.na(max.intensity)){
        next
      }
      
      # Combine data for box plot
      quantiles <- NULL ; means = NULL
      for(j in 1:length(ISI.type)){
        quantiles <- cbind(quantiles, Constant[[j]]$x[i,c(1,3,4,5,7)])
        means <- c(means,Constant[[j]]$x[i,8])
      }
      colnames(quantiles) <- c('Gam', 'IG', 'LN', 'Poi', 'Wei')
      
      pdf(paste0(new.dir,'/spike_seq_',i,'.pdf'))
      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 1, 0, 0))
      boxplot(quantiles, xlab = 'ISI type', ylab = 'Intensity', cex.lab = 2 , cex.axis = 1.5)
      points(1:length(ISI.type), means,col=2,pch=20)
      grid(nx=0,ny=10)
      dev.off()
      
    }
    
  }
  
  
}

## - Plot Model Intensity comparisons. 
{
  ISI.type <- names(PWC)
  dir.create(paste0(directory, '/Model_intensity_comparisions'), showWarnings = F)
  for(i in 1:length(ISI.type)){
    new.dir <- paste0(directory, '/Model_intensity_comparisions/',ISI.type[i])
    dir.create(new.dir, showWarnings = F)
  
      for(ii in 1:ncol(spikes)){
      # If we have spikes do:
      if(!is.na(spikes[1,ii])){
        # Work out the maximum height required in the plot. 
        max.intensity <- -10
        max.intensity <- max(max.intensity, PWC[[i]]$x_high[ii,],GP[[i]]$x_high[ii,], Constant[[i]]$x[ii,6], na.rm = T)
        
        # If we don't have any intensities break out of the current loop. 
        if(max.intensity<0 || is.na(max.intensity)){
          next
        }
        
        # Get spikes and end time.
        s <- spikes[,ii][!is.na(spikes[,ii])]
        end.time <- max(s) ; s <- s[-length(s)]
        
        colo_shade <-viridis::viridis(3,alpha=0.5)  
        colo <-viridis::viridis(3,alpha=1)  
        
        pdf(paste0(new.dir,'/spike_seq_',ii,'.pdf'))
        mar.default <- c(5,4,4,2) + 0.1
        par(mar = mar.default + c(0, 1, 0, 0))
        t <- seq(0,end.time,end.time/2000)
        plot(0,0, type='n', xlab= 'Time (s)', ylab = 'Intensity', cex.lab = 2 , cex.axis = 1.8, ylim = c(0,max.intensity), xlim = c(0,end.time))
        polygon(c(t,rev(t)), c(PWC[[i]]$x_low[ii,],rev(PWC[[i]]$x_high[ii,])), col = colo_shade[1], border=NA )
        polygon(c(t,rev(t)), c(GP[[i]]$x_low[ii,],rev(GP[[i]]$x_high[ii,])), col = colo_shade[2], border=NA )
        polygon(c(0,end.time,end.time,0), c(rep(Constant[[i]]$x[ii,2],2),rep(Constant[[i]]$x[ii,6],2)), col = colo_shade[3], border=NA )
        
        lines(t, PWC[[i]]$x_mean[ii,],col=colo[1],lwd=4)
        lines(t, GP[[i]]$x_mean[ii,],col=colo[2],lwd=4)
        lines(c(0,end.time), rep(Constant[[i]]$x[ii,8],2),col=colo[3],lwd=4)
        rug(s,ticksize=0.05,lwd=2)
        legend('topright', col = colo, legend = c('PWC', 'GP','Constant'), lty = rep(1,length(ISI.type)))
        dev.off()
        
      }
      
    }
    
  }
}

## - Model assess plots - KS/QQ + Sample sequences. 
{
  source('~/Package/R/Analysis.R')
  source('~/Package/R/QQ_KS.R')
  source('~/Package/R/Shared.R')
  filepath <- paste0(directory,'/Results.Rdata')
  KK_QQ_PWC <- get_QQ_KS(filepath, model = 'PWC', return_val = T, min.spikes= 18, ignore = ignore, rescale=T, opthyp = T)
  KK_QQ_Con <- get_QQ_KS(filepath, model = 'Constant', return_val = T, min.spikes= 18, ignore = ignore,rescale=T, opthyp=T)
  KK_QQ_GP <- get_QQ_KS(filepath, model = 'GP', return_val = T, min.spikes= 18, ignore = ignore,rescale=T, opthyp = T)
  # Gets the slopes, convert to angle.
  slopes_PWC <- KK_QQ_PWC$slopes
  slopes_Con <- KK_QQ_Con$slopes
  slopes_GP <- KK_QQ_GP$slopes
  for(i in 1:2){
    slopes_PWC[[i]] <- atan(slopes_PWC[[i]])
    slopes_Con[[i]] <- atan(slopes_Con[[i]])
    slopes_GP[[i]] <- atan(slopes_GP[[i]])
  }
  
  save(KK_QQ_PWC, KK_QQ_Con,KK_QQ_GP,slopes_GP,slopes_Con,slopes_PWC, file =paste0(directory,'/QQ_KS_results.Rdata') )
  
  # Plot simulated spikes together with spikes used to fit the model. 
  {
    source('~/Package/R/SimulatingSpikes.R')
    models <- c('Constant','PWC','GP')
    ISIs <- c("Gamma", "InverseGaussian", "LogNormal", "Exponential", "Weibull", "Exponential_Tmin")
    new.dir <- paste0(directory,'/Simulated_plots')
    dir.create(new.dir, showWarnings = F)
    for(i in 1:ncol(spikes)){
      for(kk in 1:length(models)){
        
      cat('i=',i,' and kk = ', kk,'\n')
        pdf(paste0(new.dir,'/',models[kk],'_spikeseq',i,'.pdf'))
        count=0
        for(jj in 1:6){
          print(jj)
          if(models[kk] == 'Constant'){
            int.cur <- rep(get(models[kk])[[jj]]$x[i,8],2001) 
          }
          else{
            int.cur <- get(models[kk])[[jj]]$x_mean[i,] 
          }
         
           hyp <- get(models[kk])[[jj]]$ISI_param[i,8] ; T.min.cur <- 0 ; ISI.cur <- ISIs[jj]
          end.time.real <- max(spikes[,i],na.rm=T)
          if(ISI.cur == "Exponential_Tmin"){
            ISI.cur <- "Exponential" ; T.min.cur <- hyp ; hyp = NA  ; T.min.cur <- T.min.cur*20/end.time.real
          } 
          int.cur <- int.cur*end.time.real/20
          int.cur <- inc.fineness(int.cur,10)
          set.seed(15)
          
          # Do we have any data for this sequence?
          if(is.na(int.cur[1])){next;count=count+1}
          
          simulated.spikes <- simulate_spikes(end.time=20, int.fn = int.cur, hyper=hyp, steps =2000, T.min = T.min.cur, ISI.type = ISI.cur, sequences = 15, add.end = F,do.log = T)
          simulated.spikes <- simulated.spikes * end.time.real/20
          # Plot the simulated spikes.
          
          # First create the plot
          plot(0,0,type = "n", ylim = c(1,12), xlim = c(0,max(spikes[!is.na(spikes)])),
               xlab = 'Time (s)',
               yaxt = "n", ylab = "", cex.lab = 1.8, cex.axis = 1.6, main=ISI.cur)
          
          # Add spikes the model was fitted to.
          s <- spikes[,i][!is.na(spikes[,i])] ; s <- s[-length(s)]
          points(s, rep(1,length(s)), pch = 20,cex=1.3,col=2)
          # Add the rastor of simulated spikes
          for(iii in 1:ncol(simulated.spikes)){
            spikes.cur <-  simulated.spikes[!is.na( simulated.spikes[,iii]),iii]
            y.axis <- rep(iii+1,length(spikes.cur))
            points(spikes.cur, y.axis, pch = 20,cex=1.3)
          }
        } 
        dev.off()
        if(count == 6){
          
        }
      }
     
  
      
    }
  }
  
  # Plot the ordered quantiles.
  {
    new.dir <- paste0(directory,'/Ordered_Quantile_plots')
    dir.create(new.dir, showWarnings = F)
    files <- c('KK_QQ_PWC', 'KK_QQ_Con', 'KK_QQ_GP')
    nam <-names(KK_QQ_PWC$QQ_KS)
    for(i in 1:length(files)){
      for(j in 1:ncol(get(files[i])$QQ_KS[[1]]$quantiles$exper)){
        cur <- rbind(get(files[i])$QQ_KS[[1]]$quantiles$exper[,j], get(files[i])$QQ_KS[[2]]$quantiles$exper[,j], get(files[i])$QQ_KS[[3]]$quantiles$exper[,j], get(files[i])$QQ_KS[[4]]$quantiles$exper[,j], get(files[i])$QQ_KS[[5]]$quantiles$exper[,j], get(files[i])$QQ_KS[[6]]$quantiles$exper[,j])
        pdf(paste0( new.dir,'/',files[i],'_spikes_',j,'.pdf'))
        par(mfrow=c(3,2))
        for(ii in 1:6){
          plot(cur[ii,][!is.na(cur[ii,])], xlab = 'Quantile Number', ylab = 'Exper Quantile', main = nam[ii])
        }
        dev.off()
      }
    }
  }
  
  # Violin plots
  {
    new.dir <- paste0(directory,'/QQ_KS_violin')
    dir.create(new.dir, showWarnings = F)
    
    pdf(paste0(new.dir,'/QQ_slope_PWC.pdf'))
    vioplot::vioplot(slopes_PWC$QQ)
    abline(h=pi/4, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/KS_slope_PWC.pdf'))
    vioplot::vioplot(slopes_PWC$KS)
    abline(h=pi/4, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/QQ_slope_GP.pdf'))
    vioplot::vioplot(slopes_GP$QQ)
    abline(h=pi/4, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/KS_slope_GP.pdf'))
    vioplot::vioplot(slopes_GP$KS)
    abline(h=pi/4, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/QQ_slope_Con.pdf'))
    vioplot::vioplot(slopes_Con$QQ)
    abline(h=pi/4, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/KS_slope_Con.pdf'))
    vioplot::vioplot(slopes_Con$KS)
    abline(h=pi/4, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/KS_slope.pdf'), width=15, height=5)
    slopes <- cbind(slopes_Con$KS, slopes_PWC$KS, slopes_GP$KS)
    colnames(slopes) <- c(paste0('C_',colnames(slopes_Con$KS)), paste0('P_',colnames(slopes_Con$KS)),paste0('G_',colnames(slopes_Con$KS)))
    vioplot::vioplot(slopes)
    abline(h=pi/4, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/QQ_slope.pdf'), width=15, height=5)
    slopes <- cbind(slopes_Con$QQ, slopes_PWC$QQ, slopes_GP$QQ)
    colnames(slopes) <- c(paste0('C_',colnames(slopes_Con$QQ)), paste0('P_',colnames(slopes_Con$QQ)),paste0('G_',colnames(slopes_Con$QQ)))
    vioplot::vioplot(slopes)
    abline(h=pi/4, lty='dotted')
    dev.off()
  }
  
  
# Original (boxplots)
  {
    new.dir <- paste0(directory,'/QQ_KS')
    dir.create(new.dir, showWarnings = F)
    
    pdf(paste0(new.dir,'/QQ_slope_PWC.pdf'))
    boxplot(KK_QQ_PWC$slopes$QQ)
    abline(h=1, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/KS_slope_PWC.pdf'))
    boxplot(KK_QQ_PWC$slopes$KS)
    abline(h=1, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/QQ_slope_GP.pdf'))
    boxplot(KK_QQ_GP$slopes$QQ)
    abline(h=1, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/KS_slope_GP.pdf'))
    boxplot(KK_QQ_GP$slopes$KS)
    abline(h=1, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/QQ_slope_Con.pdf'))
    boxplot(KK_QQ_Con$slopes$QQ)
    abline(h=1, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/KS_slope_Con.pdf'))
    boxplot(KK_QQ_Con$slopes$KS)
    abline(h=1, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/KS_slope.pdf'), width=15, height=5)
    slopes <- cbind(KK_QQ_Con$slopes$KS, KK_QQ_PWC$slopes$KS, KK_QQ_GP$slopes$KS)
    colnames(slopes) <- c(paste0('C_',colnames(KK_QQ_Con$slopes$KS)), paste0('P_',colnames(KK_QQ_Con$slopes$KS)),paste0('G_',colnames(KK_QQ_Con$slopes$KS)))
    boxplot(slopes, cex.axis = 0.9, ylim = c(0,4))
    abline(h=1, lty='dotted')
    dev.off()
    
    pdf(paste0(new.dir,'/QQ_slope.pdf'), width=15, height=5)
    slopes <- cbind(KK_QQ_Con$slopes$QQ, KK_QQ_PWC$slopes$QQ, KK_QQ_GP$slopes$QQ)
    colnames(slopes) <- c(paste0('C_',colnames(KK_QQ_Con$slopes$KS)), paste0('P_',colnames(KK_QQ_Con$slopes$KS)),paste0('G_',colnames(KK_QQ_Con$slopes$KS)))
    boxplot(slopes, cex.axis = 0.9, ylim = c(0,3))
    abline(h=1, lty='dotted')
    dev.off()
  }
  
}
 
 
 ### All fit into an R.file 
 directory <- getwd()
  {
 
 # 1) Collect all files into one .RData file. 
 rm(list = ls())
    nspikes  = 5
 # Store the pwc results in one variable, PWC. 

 ISI.list <- c('Gamma', 'InverseGaussian','LogNormal','Exponential','Weibull', 'Exponential_Tmin')
 PWC <- list('Gamma'=NULL, 'InverseGaussian'=NULL,'LogNormal'=NULL,'Exponential'=NULL,'Weibull'=NULL, 'Exponential_Tmin' = NULL)
 ISI <-  c('Gamma', 'NewIG','LN','Poisson','Weibull', 'PoiTmin')
 for(j in 1:length(ISI)){
   # Fill 4 tables: x_mean, x_low,x_high, ISIparam
   x_mean <- matrix(NA, ncol=2001, nrow = nspikes)
   x_low <- matrix(NA, ncol=2001, nrow = nspikes)
   x_high <- matrix(NA, ncol=2001, nrow = nspikes)
   ISI_param <- matrix(NA,ncol=8,nrow = nspikes)
   
   for(i in 1:nspikes){
     # Fill 4 tables: x_mean, x_low,x_high, ISIparam
     if(file.exists(paste0('PWC/Index',ISI[j],i,'.Rdata'))){
       load(paste0('PWC/Index',ISI[j],i,'.Rdata'))
       x_mean[i,] <- x.mean$mean
       x_low[i,] <- x.mean$lower
       x_high[i,] <- x.mean$upper
       ISI_param[i,] <- hyper.result 
       if(ISI[j] == 'PoiTmin'){
         ISI_param[i,] = tmin.result
       }
     }
   }
   cur <- list(x_mean = x_mean, x_low=x_low,x_high = x_high, ISI_param=ISI_param)
   PWC[[j]] <- cur 
 }
 
 # store the Constant results in one variable, Constant.
 ISI.list <- c('Gamma', 'InverseGaussian','LogNormal','Exponential','Weibull')
 Constant <- list('Gamma'=NULL, 'InverseGaussian'=NULL,'LogNormal'=NULL,'Exponential'=NULL,'Weibull'=NULL,'Exponential_Tmin' = NULL)
 ISI <-  c('Gamma', 'NewIG','LN','Poisson','Weibull','PoiTmin')
 for(j in 1:length(ISI)){
   # Fill 4 tables: x_mean, x_low,x_high, ISIparam
   x <- matrix(NA, ncol=8, nrow = nspikes)
   ISI_param <- matrix(NA,ncol=8,nrow = nspikes)
   
   for(i in 1:nspikes){
     # Fill 4 tables: x_mean, x_low,x_high, ISIparam
     if(file.exists(paste0('Constant/Cons',ISI[j],i,'.Rdata'))){
       load(paste0('Constant/Cons',ISI[j],i,'.Rdata'))
       x[i,] <- x.result
       ISI_param[i,] <- hyper.result 
       if(ISI[j] == 'PoiTmin'){
         ISI_param[i,] = tmin.result
       }
     }
   }
   cur <- list(x = x, ISI_param=ISI_param)
   Constant[[j]] <- cur 
 }
 
 # Store the GP results in one variable, GP. 
 ISI.list <- c('Gamma', 'InverseGaussian','LogNormal','Exponential','Weibull')
 GP <- list('Gamma'=NULL, 'InverseGaussian'=NULL,'LogNormal'=NULL,'Exponential'=NULL,'Weibull'=NULL,'Exponential_Tmin' = NULL)
 ISI <-  c('Gamma', 'NewIG','LN','Poisson','Weibull','PoiTmin')
 for(j in 1:length(ISI)){
   # Fill 4 tables: x_mean, x_low,x_high, ISIparam
   x_mean <- matrix(NA, ncol=2001, nrow = nspikes)
   x_low <- matrix(NA, ncol=2001, nrow = nspikes)
   x_high <- matrix(NA, ncol=2001, nrow = nspikes)
   ISI_param <- matrix(NA,ncol=8,nrow = nspikes)
   
   for(i in 1:nspikes){
     # Fill 4 tables: x_mean, x_low,x_high, ISIparam
     if(file.exists(paste0('GP/Index',ISI[j],i,'.Rdata'))){
       load(paste0('GP/Index',ISI[j],i,'.Rdata'))
       x_mean[i,] <- x.mean$mean  ; x_mean[i,2001] <- x_mean[i,2000]
       x_low[i,] <- x.mean$lower ; x_low[i,2001] <- x_low[i,2000]
       x_high[i,] <- x.mean$upper ; x_high[i,2001] <- x_high[i,2000]
       ISI_param[i,] <- hyper.result
       if(ISI[j] == 'PoiTmin'){
         ISI_param[i,] = tmin.result
       }
     }
   }
   cur <- list(x_mean = x_mean, x_low=x_low,x_high = x_high, ISI_param=ISI_param)
   GP[[j]] <- cur
 }
All <- list(Constant = Constant, PWC=PWC, GP=GP) 
save(All, file='results_all.RData')

 
 } 
 
