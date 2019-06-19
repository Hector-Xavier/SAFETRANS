library(plotrix)

reading_and_preparation <- function(filename,significant_data_sets)
{
  file_under_scrutiny <- file(filename, "rb")
  metadata <- rawToChar(readBin(file_under_scrutiny,raw(),n=80))
  temp_m <- rawToChar(readBin(file_under_scrutiny,raw(),n=1))
  while (is.na(strsplit(temp_m,split="\r")[[1]][2]=="\n"))
    temp_m <- paste(temp_m,rawToChar(readBin(file_under_scrutiny,raw(),n=1)),sep="")
  metadata <- c(metadata,temp_m)
  rm(temp_m)
  metadata <- c(metadata,rawToChar(readBin(file_under_scrutiny,raw(),n=80)))
  metadata <- c(metadata,rawToChar(readBin(file_under_scrutiny,raw(),n=80)))
  for (i in 5:(4+as.integer(strsplit(metadata[4],split=" ")[[1]][6])))
    metadata <- c(metadata,rawToChar(readBin(file_under_scrutiny,raw(),n=80)))
  bin_number <- as.integer(strsplit(metadata[5],split=" ")[[1]][5])
  raw_data=matrix(nrow=bin_number,ncol=significant_data_sets)
  headers <- c()
  trimming <- 0
  for (i in 1:significant_data_sets)
  {
    readBin(file_under_scrutiny,raw(),n=2)
    raw_data[,i] <- readBin(file_under_scrutiny,integer(),n=bin_number)
    headers <- c(headers,paste(c("analog","counting")[as.integer(strsplit(metadata[4+i],split=" ")[[1]][3])+1],strsplit(metadata[4+i],split=" ")[[1]][9],sep="."))
    if (as.integer(strsplit(metadata[4+i],split=" ")[[1]][12])>0)
      raw_data[,i] <- c(raw_data[(as.integer(strsplit(metadata[4+i],split=" ")[[1]][12])+1):dim(raw_data)[1],i],rep(0,as.integer(strsplit(metadata[4+i],split=" ")[[1]][12])))
    if (as.logical(as.integer(strsplit(metadata[4+i],split=" ")[[1]][3])))
    {
      raw_data[,i] <- raw_data[,i]*150/as.integer(strsplit(metadata[4+i],split=" ")[[1]][15])/as.numeric(strsplit(metadata[4+i],split=" ")[[1]][8])
    } else {
      raw_data[,i] <- raw_data[,i]*as.numeric(strsplit(metadata[4+i],split=" ")[[1]][16])*1000/as.integer(strsplit(metadata[4+i],split=" ")[[1]][15])/2^as.numeric(strsplit(metadata[4+i],split=" ")[[1]][14])
    }
    trimming <- max(trimming,as.integer(strsplit(metadata[4+i],split=" ")[[1]][12]))
  }
  close(file_under_scrutiny)
  rm(file_under_scrutiny)
  raw_data <- raw_data[1:(dim(raw_data)[1]-trimming),]
  colnames(raw_data) <- headers
  return(raw_data)
}

background_subtraction <- function(data, pre_trigger=FALSE, unsupervised=FALSE, first_bin=NULL, last_bin=NULL, first_valid_bin=NULL,bin_width=NULL,verbose=FALSE)
{
  if (unsupervised)
  {
    return(background_subtraction(data,first_bin=upper_atmosphere_cutoff(data,bin_width,verbose=verbose),last_bin=length(data)))
  } else {
    data <- data - mean(data[first_bin:last_bin])
    if (pre_trigger)
      data[1:(length(data)-first_valid_bin+1)] <- c(data[first_valid_bin:length(data)],rep(0,data[(length(data)-first_valid_bin+2):length(data)]))
    data[data<0] <- 0
  }
  return(data)
}

upper_atmosphere_cutoff <- function(data,bin_width,range_offset=0,verbose=FALSE)
{
  if (is.null(bin_width))
    stop("Please enter a valid value for bin_width, as it is necessary for unsupervised usage.")
  temp_data <- expm1(range_correction(data=background_subtraction(data,first_bin=length(data)-2000,last_bin=length(data)),bin_width))
  temp_snr <- c()
  for (i in 1:(length(temp_data)-150))
    temp_snr <- c(temp_snr,mean(temp_data[i:(i+150)])/sd(temp_data[i:(i+150)]))
  cutoff <- NULL
  for (i in 1:4000)
    if (sd(temp_snr[i:(i+600)]) < sd(temp_snr[ceiling(37500/bin_width):length(temp_snr)]))
    {
      cutoff <- i+1
      break
    }
  cutoff <- min(max(floor(8000/bin_width),cutoff),floor(20000/bin_width))
  if (verbose)
    message("Optimal cutoff height: ",cutoff*bin_width," m.")
  return(cutoff)
}

anomaly_detection <- function(data,bin_width=NULL,verbose=FALSE,confidence=0.8,window_size=20)
{
  data <- filter(data,rep(1,5)/5)
  data[is.na(data)] <- 0
  data <- as.vector(data)
  indicator <- rep(0,length(data))
  beginning <- 0
  for (i in 1:(length(data)-window_size+1-19))
  {
    if (i>1)
    {
      if (beginning==0)
        if (temp_beginning < 0.8 && sum(data[(i+1):(i+window_size-1)] > data[i:(i+window_size-2)])/window_size >= 0.8)
        {
          indicator[i] <- 1
          beginning <- i
        } 
      
      if(beginning>0)
        if (data[i] <= 3*data[beginning] && temp_ending >= 0.55 && sum(data[(i+1):(i+window_size-1)] < data[i:(i+window_size-2)])/window_size < 0.55)
        {
          indicator[i+19] <- -1
          beginning <- 0
        }
    }
    temp_beginning <- sum(data[(i+1):(i+window_size-1)] > data[i:(i+window_size-2)])/window_size
    temp_ending <- sum(data[(i+1):(i+window_size-1)] < data[i:(i+window_size-2)])/window_size
  }
  if (beginning>0)
    indicator[i-1] <- -1
  if (verbose && sum(abs(indicator))!=0)
  {
    if (is.null(bin_width))
      stop("Please enter a valid value for bin_width, as it is necessary for printing out exact heights.")
    for (i in 1:sum(indicator==1))
    {
      message("Anomaly #",i," bottom at an altitude of ",seq(length(indicator))[indicator==1][i]*bin_width," m.")
      message("Anomaly #",i," top at an altitude of ",seq(length(indicator))[indicator==-1][i]*bin_width," m.")
      message("Anomaly #",i," width estimated at ",(seq(length(indicator))[indicator==-1][i]-seq(length(indicator))[indicator==1][i])*bin_width," m.")
    }
  }
  return(indicator)
}

range_correction <- function(data,bin_width,range_offset=0)
{
  return(log1p(data*(seq(length(data))*bin_width+range_offset)))
}

extinction_coefficient <- function(data,bin_width,k=1,verbose=FALSE)
{
  corrected_data <- range_correction(background_subtraction(data=data,unsupervised=TRUE,bin_width=bin_width,verbose=verbose),bin_width=bin_width)
  anomalies <- anomaly_detection(data=expm1(corrected_data),bin_width=bin_width,verbose=verbose)
  if(sum(abs(anomalies))==0)
    anomalies[length(anomalies)] <- 1
  cutoff <- max(min(seq(length(data))[anomalies==1][1],upper_atmosphere_cutoff(data,bin_width=bin_width)),600)
  indices <- ((cutoff-199):cutoff)[corrected_data[(cutoff-199):cutoff]!=0]
  signal <- corrected_data[indices]
  coefficient_estimate <- (signal[1:(length(signal)-1)]-signal[2:length(signal)])/(indices[2:length(indices)]-indices[1:(length(indices)-1)])/bin_width/2
  coefficient_estimate <- mean(coefficient_estimate[coefficient_estimate>0])
  extinction <- exp((corrected_data[min(seq(length(anomalies))[anomalies==1][1],seq(50)[corrected_data[1:50]==max(corrected_data[1:50])]):length(corrected_data)]-mean(signal))/k)
  extinction[extinction<0] <- 0
  integrals <- c()
  for (i in 1:length(extinction))
    integrals <- c(integrals,sum(extinction[i:cutoff])*bin_width)
  return(c(rep(0,min(seq(length(anomalies))[anomalies==1][1],seq(50)[corrected_data[1:50]==max(corrected_data[1:50])])-1),extinction/(coefficient_estimate^(-1)+2*integrals/k)))
}

scanning_profile_extinction <- function(scanning_directory,measurements_of_interest,is_scan,k=1,verbose=FALSE,output_file=FALSE)
{
  snr_limit <- 0.05
  file_list <- list.files(scanning_directory)[list.files(scanning_directory)!="Output"]
  temp <- reading_and_preparation(filename=paste(scanning_directory,file_list[1],sep="/"),significant_data_sets=max(measurements_of_interest,2))[,measurements_of_interest]
  if (length(measurements_of_interest)==1)
  {
    data <- matrix(ncol=1,nrow=length(temp))
    data[,1] <- temp
  } else {
    data <- matrix(ncol=length(measurements_of_interest),nrow=dim(temp)[1])
    data[,1:length(measurements_of_interest)] <- temp
  }
  bin_width <- as.numeric(read.table(paste(scanning_directory,file_list[1],sep="/"),skip=4,nrows=1,fill=TRUE)[7])
  for (i in 2:length(file_list))
    data <- cbind(data,reading_and_preparation(filename=paste(scanning_directory,file_list[i],sep="/"),significant_data_sets=max(measurements_of_interest,2))[,measurements_of_interest])
  if (!is_scan)
  {
    message("Processing data as separate measurements.")
  } else {
    for (i in 1:length(measurements_of_interest))
    {
      indices_under_scrutiny <- seq(dim(data)[2])[(seq(dim(data)[2])-i)%%length(measurements_of_interest)==0]
      if (sd(data[1,indices_under_scrutiny]*bin_width^2)/mean(data[1,indices_under_scrutiny]*bin_width^2)>snr_limit)
      {
        message("Measurements in channel ",measurements_of_interest[i]," at a distance of ",bin_width," m exhibit SNR less than ",ceiling(1/snr_limit),".")
        message("Current SNR of channel ",measurements_of_interest[i],": ",floor(mean(data[1,indices_under_scrutiny]*bin_width^2)/sd(data[1,indices_under_scrutiny]*bin_width^2)),". Caution and possible instrument calibrartion recommended.")
      }
      stabilization <- c()
      for (j in 1:length(indices_under_scrutiny))
        stabilization <- c(stabilization,mean(data[1,indices_under_scrutiny])/data[1,indices_under_scrutiny[j]])
      data[,indices_under_scrutiny] <- t(t(data[,indices_under_scrutiny])*stabilization)
    }    
  }
  if (!verbose)
    progress <- txtProgressBar(max=dim(data)[2],char="=",style=3)
  for (i in 1:dim(data)[2])
  {
    if (verbose && (i%%length(measurements_of_interest)==1 || length(measurements_of_interest)==1))
      message(ceiling(i/length(measurements_of_interest)),"/",length(file_list)," : Processing data file ",file_list[ceiling(i/length(measurements_of_interest))])
    data[,i] <- extinction_coefficient(data=data[,i],bin_width,k,verbose)
    if (!verbose)
      setTxtProgressBar(progress,i)
  }
  if (!verbose)
    close(progress)
  headers <- c()
  for (i in 1:length(file_list))
    headers <- c(headers,paste(rep(paste("Zenith",-as.numeric(read.table(paste(scanning_directory,file_list[i],sep="/"),skip=1,nrows=1,fill=TRUE)[9]),"Azimuth",(as.numeric(read.table(paste(scanning_directory,file_list[i],sep="/"),skip=1,nrows=1,fill=TRUE)[10])+as.numeric(read.table(paste(scanning_directory,file_list[i],sep="/"),skip=2,nrows=1,fill=TRUE))[7])%%360,sep="_"),length(measurements_of_interest)),measurements_of_interest,sep="_"))
  colnames(data) <- headers
  rownames(data) <- as.character(seq(dim(data)[1])*bin_width)
  if(output_file)
    write.table(data,file=paste("Radial_extinction_coefficients_1054_nm.txt",sep="/"),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
  return(data)
}

visibility_range <- function(extinction,bin_width,model=NULL,incoming=FALSE,incoming_range=NULL,verbose=FALSE)
{
  if (!is.null(model))
    if (sum(model==c("urban-rural","maritime"))==0)
      stop("Incorrect model designation. Correct designations are: unrban-rural, maritime, NULL")
  if (incoming && is.null(incoming_range))
    stop("Please provide the range of the incoming object in metres.")
  if (verbose)
  {
    if (!is.null(model))
    {
      message("Selected model: ",as.character(model))
      message("Translating visibility from 1054 nm to visible range.")
    } else {
      message("No model selected. Calculating visibility for a wavelength of 1054 nm.")
    }
  }
  if (!is.null(model))
  {
    if (model=="urban-rural")
      extinction <- (extinction/0.314)^(1/1.11)
    if (model=="maritime")
      extinction <- (extinction/0.996)^(1/1.20)
  }
  visibility <- c()
  if (incoming)
  {
    for (i in ceiling(incoming_range/bin_width):1)
      visibility <- c(visibility,sum(extinction[ceiling(incoming_range/bin_width):i]*bin_width))
    if (visibility[length(visibility)]>=3)
    {
      return(length(visibility[visibility<=3])*bin_width)
    } else {
      return(incoming_range)
    }
  } else {
    for (i in 1:length(extinction))
      visibility <- c(visibility,sum(extinction[1:i]*bin_width))
    return(length(visibility[visibility<=3])*bin_width)
  }
}

radial_visibility_profile <- function(extinction_profile,is_scan=TRUE,model=NULL,output_file=TRUE,verbose=FALSE)
{
  bin_width <- as.numeric(rownames(extinction_profile)[2])-as.numeric(rownames(extinction_profile)[1])
  visibility <- c()
  if (verbose)
  {
    for (i in 1:dim(extinction_profile)[2])
    {
      message(i,"/",dim(extinction_profile)[2]," - Measurement set designation: ",colnames(extinction_profile)[i])
      visibility <- c(visibility,visibility_range(extinction=extinction_profile[,i],bin_width,model=model,verbose=TRUE))
      message("Outward visibility: ",visibility[i]," m.")
    }
  } else {
    progress <- txtProgressBar(max=dim(extinction_profile)[2],char="=",style=3)
    for (i in 1:dim(extinction_profile)[2])
    {
      visibility <- c(visibility,visibility_range(extinction=extinction_profile[,i],bin_width,model=model))
      setTxtProgressBar(progress,i)
    }
    close(progress)
  }
  if (is.null(model))
    model <- "no_model"
  if (output_file)
    if(!file.exists(paste("Radial_outward_visibility_distance_",model,".txt",sep="")))
      write.table(visibility,file=paste("Radial_outward_visibility_distance_",model,".txt",sep=""),quote=FALSE,sep="\t",row.names=colnames(extinction_profile),col.names="Visibility_in_metres")
  visibility <- matrix(visibility,ncol=1,dimnames=list(colnames(extinction_profile),model))
  if (output_file && is_scan)
  {
    angle <- c()
    channels <- c()
    for (i in 1:length(visibility))
    {
      angle <- c(angle,as.integer(strsplit(rownames(visibility)[i],split="_")[[1]][4]))
      channels <- c(channels,as.integer(strsplit(rownames(visibility)[i],split="_")[[1]][5]))
    }
    
    if (length(levels(as.factor(channels)))>1)
    {
      for (i in 1:length(levels(as.factor(channels))))
      {
        png(file=file.path(getwd(),paste("Radial_visibility_",model,"_channel_",levels(as.factor(channels))[i],".png", sep = "")),width=1000,height=1000)
        radial.plot(lengths = visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1],radial.pos = angle[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i]/180*pi,labels=c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW"),label.pos=((seq(16)-1)*22.5-angle[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i][1])/180*pi,start=-(angle[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i][1]/180-1/2)*pi,clockwise=TRUE,rp.type="p",radial.lim = pretty(range(visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1])),show.grid.labels=1,radial.labels = paste(pretty(range(visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1]))/1000,"km",sep=" "),show.centroid=TRUE,main=paste("Radial visibility at",strsplit(rownames(visibility)[1],split="_")[[1]][2],"degrees in channel",levels(as.factor(channels))[i],sep=" "))
        dev.off()
      }
    } else {
      png(file=file.path(getwd(),paste("Radial_visibility_",model,".png", sep = "")),width=1000,height=1000)
      radial.plot(lengths = visibility[,1],radial.pos = angle/180*pi,labels=c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW"),label.pos=((seq(16)-1)*22.5-angle[1])/180*pi,start=-(angle[1]/180-1/2)*pi,clockwise=TRUE,rp.type="p",radial.lim = pretty(range(visibility[,1])),show.grid.labels=1,radial.labels = paste(pretty(range(visibility[,1]))/1000,"km",sep=" "),show.centroid=TRUE,main=paste("Radial visibility at ",strsplit(rownames(visibility)[1],split="_")[[1]][2]," degrees [model: ",model,"]",sep=""))
      dev.off()
    }
    if(!file.exists(paste(getwd(),"Azimuth_visibility_plots",sep="/")))
      dir.create(paste(getwd(),"Azimuth_visibility_plots",sep="/"))
    if (output_file)
      message("Calculating radial visibility over visible range. Selected atmospheric model: ",model)     
    if (!verbose)
      progress <- txtProgressBar(max=length(visibility),char="=",style=3)
    for (i in 1:length(visibility))
    {
      if (verbose)
      {
        message(i,"/",dim(extinction_profile)[2]," - Calculating visibilities of measurement set: ",colnames(extinction_profile)[i])
      } else {
        setTxtProgressBar(progress,i)
      }
      if (model=="no_model")
        model <- NULL
      incoming_vis <- c()
      for (j in seq(6000,1,-50))
        incoming_vis <- c(incoming_vis,visibility_range(extinction_profile[,i],bin_width,model,incoming=TRUE,incoming_range=j*bin_width,verbose=FALSE))
      if (is.null(model))
        model <- "no_model"
      png(file=file.path(paste(paste(getwd(),"Azimuth_visibility_plots",sep="/"),paste(model,"_Visibility_",i,"_",rownames(visibility)[i],".png", sep = ""),sep="/")),width=2000,height=1600)
      barplot(incoming_vis[length(incoming_vis):1]/1000,log="y",col=c("red","grey")[(seq(length(incoming_vis))*50*bin_width>=visibility[i,1])+1],names.arg=paste(seq(6000,1,-50)[length(incoming_vis):1]*bin_width/1000,"km",sep=" "),las=2,xlab="Distance",ylab="Visibility (km)",main=paste(rownames(visibility)[i], " [model: ",model,"]",sep=""))
      legend("topleft",legend=c("Optical contact with overhead airspace","No optical contact"),pch=15,col=c("red","grey"))
      dev.off()
    }
    if (!verbose)
      close(progress)
  }
  return(visibility)
}

cartesian_visibility_profile <- function(extinction_profile,model=NULL,incoming=FALSE,incoming_distance=NULL,incoming_height=NULL,output_files=TRUE,verbose=TRUE)
{
  if (incoming && (is.null(incoming_distance) || is.null(incoming_height)))
    stop("Please designate both the height and the distance of the incoming object.")
  maximum_height <- 15000
  bin_width <- as.numeric(rownames(extinction_profile)[2])-as.numeric(rownames(extinction_profile)[1])
  cartesian_profile <- matrix(0,nrow=ceiling(maximum_height/bin_width),ncol=dim(extinction_profile)[1])
  colnames(cartesian_profile) <- paste("Distance",rownames(extinction_profile)[1:dim(cartesian_profile)[2]],"m",sep="_")
  rownames(cartesian_profile) <- paste("Height",rownames(extinction_profile)[1:dim(cartesian_profile)[1]],"m",sep="_")
  if (incoming && (ceiling(incoming_distance/bin_width)>dim(cartesian_profile)[2] || incoming_height>maximum_height))
    stop("Designated coordinates are beyond the detectable range. Maximum height: ",maximum_height," m, maximum distance: ",dim(cartesian_profile)[2]*bin_width," m.")
  angles <- c()
  for (i in 1:dim(extinction_profile)[2])
    angles <- c(angles,as.numeric(strsplit(colnames(extinction_profile)[i],split="_")[[1]][2]))
  for (i in 1:length(angles))
    for (j in 1:dim(cartesian_profile)[1])
    {
      y <- ceiling(sinpi((angles[i])/180)*seq(dim(cartesian_profile)[1]))
      if (sum(y)==0)
        y <- y+1
      x <- ceiling(cospi((angles[i])/180)*seq(dim(cartesian_profile)[2]))
      if (sum(x)==0)
        x <- x+1
      cartesian_profile[y[j],x[j]] <- extinction_profile[j,i]
    }
  progress <- txtProgressBar(max=dim(cartesian_profile)[1],char="=",style=3)
  for (i in 2:dim(cartesian_profile)[1])
  {
    if (sum(cartesian_profile[i,]>0)==1)
    {
      cartesian_profile[i,] <- cartesian_profile[i,cartesian_profile[i,]!=0]
    } else if (sum(cartesian_profile[i,]!=0) > 1) {
      indices <- seq(dim(cartesian_profile)[2])[cartesian_profile[i,]!=0]
      cartesian_profile[i,1:indices[1]] <- cartesian_profile[i,indices[1]]
      cartesian_profile[i,indices[length(indices)]:dim(cartesian_profile)[2]] <- cartesian_profile[i,indices[length(indices)]]
      for (j in 2:length(indices))
        if ((indices[j-1]+1)!=indices[j])
          cartesian_profile[i,(indices[j-1]+1):(indices[j]-1)] <- cartesian_profile[i,indices[j-1]]+(cartesian_profile[i,indices[j]]-cartesian_profile[i,indices[j-1]])*seq((indices[j]-indices[j-1]-1))/(indices[j]-indices[j-1]-1)
    }
    setTxtProgressBar(progress,i)
  }
  close(progress)
  non_zero <- seq(dim(cartesian_profile)[1])[rowSums(cartesian_profile)==0][1]-1
  hold_the_door <- c()
  for (i in seq(dim(cartesian_profile)[1])[rowSums(cartesian_profile)==0][1]:(seq(dim(cartesian_profile)[1])[rowSums(cartesian_profile)==0][sum(rowSums(cartesian_profile)==0)]+1))
  {
    if (sum(cartesian_profile[i,])==0)
    {
      if (non_zero>0)
      {
        cartesian_profile[i,] <- cartesian_profile[non_zero,]
      } else {
        hold_the_door <- c(hold_the_door,i)
      }
      
    } else {
      non_zero <- i
      if (length(hold_the_door)>0)
        for (j in 1:length(hold_the_door))
          cartesian_profile[hold_the_door[j],] <- cartesian_profile[i,]
        hold_the_door <- c()
    }  
  }
  if (!is.null(incoming_height) && !is.null(incoming_distance))
  {
    if(incoming)
    {
      horizontal_visibility <- visibility_range(extinction=cartesian_profile[ceiling(incoming_height/bin_width),1:ceiling(incoming_distance/bin_width)],bin_width,model,incoming,incoming_distance,verbose)
    } else {
      horizontal_visibility <- visibility_range(extinction=c(cartesian_profile[ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width):dim(cartesian_profile)[2]],rep(cartesian_profile[ceiling(incoming_height/bin_width),dim(cartesian_profile)[2]],ceiling(20000/bin_width))),bin_width,model,incoming,incoming_distance,verbose)
    }
    vertical_visibility <- visibility_range(extinction=cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)],bin_width,model,incoming=TRUE,incoming_height,verbose=FALSE)
    if (vertical_visibility < incoming_height)
    {
      slant_visibility <- "No optical contact between object and ground. Slant visibility unavailable."
    } else {
      pseudo_visibility <- visibility_range(extinction=c(rep(cartesian_profile[1,ceiling(incoming_distance/bin_width)],maximum_height/bin_width),cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),bin_width,model,incoming=TRUE,incoming_height+maximum_height,verbose=FALSE)
      slant_visibility <- ceiling(sqrt(pseudo_visibility^2 - incoming_height^2))
    }
    if (verbose)
    {      
      if (incoming)
      {
        if (horizontal_visibility < incoming_distance)
        {
          message("Horizontal incoming visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m : ",horizontal_visibility," m.")
        } else {
          message("Incoming object at a height of ",incoming_height," m and distance of ",incoming_distance," m has optical contact with overhead airspace.")
        }
      } else {
        message("Horizontal visibility of outcoming object at a height of ",incoming_height," m and distance of ",incoming_distance," m : ",horizontal_visibility," m.")
      }
      if (vertical_visibility < incoming_height)
      {
        message("Vertical visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m : ",vertical_visibility," m.")
        message(slant_visibility)
      } else {
        message(c("Outcoming","Incoming")[as.integer(incoming)+1]," object at a height of ",incoming_height," m and distance of ",incoming_distance," m has optical contact with ground.")
        message("Slant visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m : ",slant_visibility," m.")
      }
    }
  }
  if (output_files)
  {
    if (is.null(model))
      model <- "no_model"
    message("Writing output files to disk...")
    if (!file.exists(paste("Cartesian_extinction_profile_",model,".txt.gz",sep="")))
    {
      write.table(cartesian_profile,file=paste("Cartesian_extinction_profile_",model,".txt",sep=""),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
      system(paste("gzip ","Cartesian_extinction_profile_",model,".txt",sep=""))      
    }
    if (!is.null(incoming_height) && !is.null(incoming_distance))
    {
      if (incoming)
      {
        write.table(c(horizontal_visibility,vertical_visibility,slant_visibility),file=paste("Incoming_object_visibility_",model,".txt",sep=""),quote=FALSE,sep="\t",col.names=paste("Incoming height: ",incoming_height," m, incoming distance: ",incoming_distance," m.",sep=""),row.names=c("Horizontal visibility: ","Vertical visibility: ","Slant visibility:"))
      } else {
        write.table(c(horizontal_visibility,vertical_visibility,slant_visibility),file=paste("Outcoming_object_visibility_",model,".txt",sep=""),quote=FALSE,sep="\t",col.names=paste("Outcoming height: ",incoming_height," m, outcoming distance: ",incoming_distance," m.",sep=""),row.names=c("Horizontal visibility: ","Vertical visibility: ","Slant visibility:"))
      }
    }
  }
  #return(cartesian_profile)
}
 
