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
  
  #in case the "extra line" with scan metadata is not used, the next line is replaced with the one after it
  metadata <- c(metadata,rawToChar(readBin(file_under_scrutiny,raw(),n=80))) 
  #metadata <- c(metadata,metadata[3])
  
  metadata <- c(metadata,rawToChar(readBin(file_under_scrutiny,raw(),n=80)))
  for (i in 5:(4+as.integer(strsplit(metadata[4],split=" ")[[1]][6])))
    metadata <- c(metadata,rawToChar(readBin(file_under_scrutiny,raw(),n=80)))
  
  #the commented line may come in handy in case of file format issues
  bin_number <- as.integer(strsplit(metadata[5],split=" ")[[1]][5])
  #bin_number <- as.integer(strsplit(metadata[6],split=" ")[[1]][5])
  
  #the commented line may come in handy in case of file format issues
  wavelength <- as.integer(strsplit(strsplit(metadata[5],split=" ")[[1]][9],split=".p")[[1]][1])
  #wavelength <- as.integer(strsplit(strsplit(metadata[6],split=" ")[[1]][9],split=".p")[[1]][1])
  raw_data=matrix(nrow=bin_number,ncol=significant_data_sets)
  
  #the metadata are used in order to name the table columns and guide value conversion as per the Raymetrics manual
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

background_subtraction <- function(data, pre_trigger=FALSE,unsupervised=FALSE,first_bin=NULL,last_bin=NULL,first_valid_bin=NULL,bin_width=NULL,zeroing=FALSE,verbose=FALSE)
{
  if (unsupervised)
  {
    #unsupervised execution will use the "upper_atmosphere_cutoff" method to choose a cutoff point for background estimation and use it to re-run itself
    return(background_subtraction(data,first_bin=upper_atmosphere_cutoff(data,bin_width=bin_width,verbose=verbose),last_bin=length(data),zeroing=zeroing))
  } else {
    data <- data - mean(data[first_bin:last_bin])
    if (pre_trigger)
      data[1:(length(data)-first_valid_bin+1)] <- c(data[first_valid_bin:length(data)],rep(0,data[(length(data)-first_valid_bin+2):length(data)]))
    
    #bad solution
    if (!zeroing)
      for (i in 2: length(data))
        if (data[i] <= 0)
          data[i] <- data[i-1]
    
    #negative signal is set to zero
    data[data<0] <- 0
  }
  return(data)
}

upper_atmosphere_cutoff <- function(data,bin_width,range_offset=0,verbose=FALSE)
{
  if (is.null(bin_width))
    stop("Please enter a valid value for bin_width, as it is necessary for unsupervised usage.")
  
  #usage of the distribution's tail as a reference for average noise profile
  temp_data <- expm1(range_correction(data=background_subtraction(data,first_bin=length(data)-2000,last_bin=length(data),zeroing=TRUE),bin_width))
  temp_snr <- c()
  for (i in 1:(length(temp_data)-150))
    temp_snr <- c(temp_snr,mean(temp_data[i:(i+150)])/sd(temp_data[i:(i+150)]))
  if (sum(is.na(temp_snr))>0)
    temp_snr[is.na(temp_snr)] <- 0
  cutoff <- NULL
  
  #a metric of the standard deviation's local stability with reference to the distribution's tail is used to estimate a candidate cutoff point
  for (i in 1:(30000/bin_width))
    if (sd(temp_snr[i:(i+600)]) < sd(temp_snr[ceiling(37500/bin_width):length(temp_snr)]))
    {
      cutoff <- i+1
      break
    }
  
  #extreme cutoff distances are not acceptable
  cutoff <- min(max(floor(2000/bin_width),cutoff),floor(20000/bin_width))
  if (verbose)
    message("Optimal cutoff range: ",cutoff*bin_width," m.")
  return(cutoff)
}

anomaly_detection <- function(data,bin_width=NULL,angle=NULL,verbose=FALSE,confidence=0.8,window_size=20)
{
  #the method checks for significant signal anomalies which could signify local atmospheric phenomena, influencing the acceptable cutoff point choice
  data <- filter(data,rep(1,5)/5)
  data[is.na(data)] <- 0
  data <- as.vector(data)
  indicator <- rep(0,length(data))
  beginning <- 0
  for (i in 1:(length(data)-window_size+1-19))
  {
    #a variant of the slope method is used on the signal distribution in order to identify candidate beginnings and ends of atmospheric anomalies
    #the sliding window size and anomaly beginning/end confidences were empirically set after result comparisons with researcher's assessements of anomaly presence
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
  
  #making sure that beginnings and ends follow a logical order
  if (beginning>0)
    indicator[i-1] <- -1
  
  #printing out anomalies may come in handy in some situations
  if (verbose && sum(abs(indicator))!=0)
  {
    if (is.null(bin_width))
      stop("Please enter a valid value for bin_width, as it is necessary for printing out exact heights.")
    for (i in 1:sum(indicator==1))
    {
      message("Anomaly #",i," bottom at an altitude of ",sinpi(angle/180) * seq(length(indicator))[indicator==1][i]*bin_width," m.")
      message("Anomaly #",i," top at an altitude of ",sinpi(angle/180) * seq(length(indicator))[indicator==-1][i]*bin_width," m.")
      message("Anomaly #",i," width estimated at ",sinpi(angle/180) * (seq(length(indicator))[indicator==-1][i]-seq(length(indicator))[indicator==1][i])*bin_width," m.")
    }
  }
  return(indicator)
}

cloud_detection <- function(filename,data,angle,bin_width,verbose,output_files)
{
  #a comprehensive cloud detection algorithm, based on Rania Soupiona's (NTUA) thesis, used for comprehensive assessment of cloud coverage for field applications
  
  if (verbose)
    message("Performing cloud detection at ",filename,".")
  
  clouds <- list(bottom=NULL,top=NULL)
  if (angle == 0)
    angle <- 1
  
  #the next three commented lines (replacing the three lines followin them) could provide an automated estimation, diverging from the original algorithm
  
  #unsupervised_cutoff <- upper_atmosphere_cutoff(data,bin_width=bin_width,verbose=verbose)
  #noise <- 3 * sd(uncorrected_signal[unsupervised_cutoff:length(uncorrected_signal)])
  #mean_background <- mean(uncorrected_signal[unsupervised_cutoff:length(uncorrected_signal)])
  
  uncorrected_signal <- data[1:min(floor(18000/sin(angle/180*pi)/bin_width),length(data))]
  noise <- 3 * sd(uncorrected_signal[(length(uncorrected_signal)-500):length(uncorrected_signal)])
  mean_background <- mean(uncorrected_signal[(length(uncorrected_signal)-500):length(uncorrected_signal)])
  
  ##smoothing from Zhou et al., ommited following a number of tests. Included as a comment for possible future use and completeness
  #temporary_signal <- uncorrected_signal
  #signal_temporary <- uncorrected_signal
  #dub_step <- temporary_signal[1]
  #step_dub <- signal_temporary[length(signal_temporary)]
  #for (i in 2:length(uncorrected_signal))
  #{
  #  if (abs(temporary_signal[i] - temporary_signal[i-1]) < noise)
  #  {
  #    temporary_signal[i] <- dub_step
  #  } else {
  #    dub_step <- temporary_signal[i]
  #  }
  #  
  #  if (abs(signal_temporary[length(signal_temporary)-i+1] - signal_temporary[length(signal_temporary)-i+2]) < noise)
  #  {
  #    signal_temporary[length(signal_temporary)-i+1] <- step_dub
  #  } else {
  #    step_dub <- signal_temporary[length(signal_temporary)-i+1]
  #  }
  #}
  ##implementation of smoothing
  #uncorrected_signal[1] <- temporary_signal[1]
  #uncorrected_signal[length(uncorrected_signal)] <- signal_temporary[length(signal_temporary)]
  #uncorrected_signal[2:length(uncorrected_signal)] <- rowMeans(cbind(temporary_signal[2:length(temporary_signal)],signal_temporary[1:(length(signal_temporary)-1)]))
  #rm(temporary_signal,signal_temporary,dub_step,step_dub,i)
  
  #the following section performs signal discretization and enhancement
  sorted_signal <- sort(uncorrected_signal,decreasing=FALSE,index.return=TRUE)
  normalization_factor <- seq(from=1,to=length(sorted_signal$x))/length(sorted_signal$x)
  for (i in 2:length(normalization_factor))
  {
    if (diff(sorted_signal$x)[i-1] == 0)
      normalization_factor[i] <- normalization_factor[i-1]
  }
  enhanced_signal <- (normalization_factor * diff(range(uncorrected_signal)) + min(uncorrected_signal))[sort(sorted_signal$ix,decreasing=FALSE,index.return=TRUE)$ix]
  
  #changepoints are estimated as divergences from a linear distribution
  changepoints <- diff((enhanced_signal - (seq(from=range(uncorrected_signal)[2],to=range(uncorrected_signal)[1],length.out=length(uncorrected_signal)))) > 0)
  if (sum(changepoints == 1) > 0 && sum(changepoints == -1) > 0)
  {
    threshold <- ceiling(45/sin(angle/180*pi)/bin_width)
    #layers below ovelap and above 12km
    changepoints[1:((seq(50)[uncorrected_signal[1:50]==max(uncorrected_signal[1:50])])[1] +1 )] <- 0
    changepoints[1:(range(seq(length(changepoints))[changepoints == 1])[1]-1)] <- 0
    changepoints[min(ceiling(12000/sin(angle/180*pi)/bin_width),length(changepoints)):length(changepoints)] <- 0
    if (sum(changepoints == 1) > 0 && sum(changepoints == -1) > 0)
    {
      changepoints[(range(seq(length(changepoints))[changepoints == -1])[2]+1):length(changepoints)] <- 0
    } else {
      changepoints <- rep(0,length(changepoints))
    }
    
    #making sure changepoints follow a logical order
    in_clouds <- TRUE
    if (sum(changepoints == 1) > 0 & sum(changepoints == -1) > 0)
    {
      for (i in (range(seq(length(changepoints))[changepoints == 1])[1]+1):(range(seq(length(changepoints))[changepoints == -1])[2]-1))
      {
        if (changepoints[i] == as.numeric(in_clouds)*2-1)
          changepoints[i] <- 0
        if (changepoints[i] == as.numeric(!in_clouds)*2-1)
          in_clouds <- !in_clouds
      }
      changepoints[(range(seq(length(changepoints))[changepoints == -1])[2]+1):length(changepoints)] <- 0
    }
    #message(sum(changepoints==1))
    
    #rejecting potential clouds under a certain thickness 
    temp_changepoints <- changepoints
    if (sum(changepoints == 1) > 0 & sum(changepoints == -1) > 0)
    {
      #debug file
      #write.table(cbind(seq(length(changepoints))[changepoints == 1]*bin_width*sin(angle/180*pi),seq(length(changepoints))[changepoints == -1]*bin_width*sin(angle/180*pi)),file=paste("Cloud layers [",filename,"] - before thin rejection.txt",sep=""),quote=FALSE,sep="\t",row.names=paste("Cloud layer",seq(length(seq(length(changepoints))[changepoints == 1])),sep=" #"),col.names=c("Cloud layer bottom [height in m.]","Cloud layer top [height in m.]"))
      
      for (i in 1:length(seq(length(changepoints))[changepoints == 1]))
      {
        if ((seq(length(changepoints))[changepoints == -1][i] - seq(length(changepoints))[changepoints == 1][i]) < threshold)
        {
          #message("Thin cloud layer erased")
          temp_changepoints[seq(length(changepoints))[changepoints == -1][i]] <- 0
          temp_changepoints[seq(length(changepoints))[changepoints == 1][i]] <- 0
        }
      }
      changepoints <- temp_changepoints
    }
    
    #merging potential clouds that are too close to one another
    if (sum(changepoints == 1) > 1 & sum(changepoints == -1) > 1)
    {
      #debug file
      #write.table(cbind(seq(length(changepoints))[changepoints == 1]*bin_width*sin(angle/180*pi),seq(length(changepoints))[changepoints == -1]*bin_width*sin(angle/180*pi)),file=paste("Cloud layers [",filename,"] - before layer merging.txt",sep=""),quote=FALSE,sep="\t",row.names=paste("Cloud layer",seq(length(seq(length(changepoints))[changepoints == 1])),sep=" #"),col.names=c("Cloud layer bottom [height in m.]","Cloud layer top [height in m.]"))
      
      for (i in 2:length(seq(length(changepoints))[changepoints == -1]))
      {
        if ((seq(length(changepoints))[changepoints == 1][i] - seq(length(changepoints))[changepoints == -1][i-1]) < threshold)
        {
          #message("Merged cloud layers")
          temp_changepoints[seq(length(changepoints))[changepoints == -1][i]] <- 0
          temp_changepoints[seq(length(changepoints))[changepoints == 1][i]] <- 0
        }
      }
      changepoints <- temp_changepoints
    }
    
    #rejecting potential clouds that are composed of very low signal intensities and could potentially be artifacts
    if (sum(changepoints == 1) > 0 & sum(changepoints == -1) > 0)
    {
      #debug file
      #write.table(cbind(seq(length(changepoints))[changepoints == 1]*bin_width*sin(angle/180*pi),seq(length(changepoints))[changepoints == -1]*bin_width*sin(angle/180*pi)),file=paste("Cloud layers [",filename,"] - before noise-based rejection.txt",sep=""),quote=FALSE,sep="\t",row.names=paste("Cloud layer",seq(length(seq(length(changepoints))[changepoints == 1])),sep=" #"),col.names=c("Cloud layer bottom [height in m.]","Cloud layer top [height in m.]"))
      
      cloud_max_signal <- c()
      for (i in 1:length(seq(length(changepoints))[changepoints == 1]))
      {
        #max(uncorrected_signal[seq(length(changepoints))[changepoints == 1][i]:seq(length(changepoints))[changepoints == -1][i]]) < mean_background + 1/3*noise
        #round(max(uncorrected_signal[seq(length(changepoints))[changepoints == 1][i]:seq(length(changepoints))[changepoints == -1][i]]),3) <= round(mean_background + 1/3*noise,3)
        if (round(max(uncorrected_signal[seq(length(changepoints))[changepoints == 1][i]:seq(length(changepoints))[changepoints == -1][i]]),3) <= round(mean_background + 1/3*noise,3))
        {
          #message("Removed low-signal cloud layer")
          #message(max(uncorrected_signal[seq(length(changepoints))[changepoints == 1][i]:seq(length(changepoints))[changepoints == -1][i]])," versus a threshold of: ",mean_background + 2/3*noise)
          
          temp_changepoints[seq(length(changepoints))[changepoints == -1][i]] <- 0
          temp_changepoints[seq(length(changepoints))[changepoints == 1][i]] <- 0
        } else {
          #message(max(uncorrected_signal[seq(length(changepoints))[changepoints == 1][i]:seq(length(changepoints))[changepoints == -1][i]])," versus a threshold of: ",mean_background + 2/3*noise)
          cloud_max_signal <-c(cloud_max_signal,max(uncorrected_signal[seq(length(changepoints))[changepoints == 1][i]:seq(length(changepoints))[changepoints == -1][i]]))
        }
      }
      changepoints <- temp_changepoints
    }
    rm(i,in_clouds,temp_changepoints)
    
    #list of cloud coverage as a vector with "1" signifying the bottom of a cloud and "-1" the top
    clouds <- list(bottom=seq(length(changepoints))[changepoints == 1],top=seq(length(changepoints))[changepoints == -1])
  }
  
  #it may be useful to see a list of identified clouds as a realtime message referencing their height and thickness
  if (verbose && length(clouds$bottom) > 0 && length(clouds$bottom) == length(clouds$top))
    for (i in 1:length(clouds$bottom))
      message("Cloud layer #",i,", from ",clouds$bottom[i]*bin_width*sin(angle/180*pi)," m to ",clouds$top[i]*bin_width*sin(angle/180*pi)," m and a cloud depth of ",(clouds$top[i]-clouds$bottom[i])*bin_width*sin(angle/180*pi),"m.")
  
  if (output_files && length(clouds$bottom) > 0 && length(clouds$bottom) == length(clouds$top))
  {
    if(!file.exists(paste(getwd(),"Cloud_layers",sep="/")))
      dir.create(paste(getwd(),"Cloud_layers",sep="/"))
    write.table(cbind(clouds$bottom*bin_width*sin(angle/180*pi),clouds$top*bin_width*sin(angle/180*pi),cloud_max_signal,rep(mean_background + 1/3*noise,length(cloud_max_signal))),file=paste(getwd(),"Cloud_layers",paste("Cloud layers [",filename,"].txt",sep=""),sep="/"),quote=FALSE,sep="\t",row.names=paste("Cloud layer",seq(length(clouds$bottom)),sep=" #"),col.names=c("Cloud layer bottom [height in m.]","Cloud layer top [height in m.]","Cloud maximum signal","Signal cutoff for noise"))
  }
  #this return statement was used only during debugging, the code normally writes a file with the results
  #return(clouds)
}

range_correction <- function(data,bin_width,range_offset=0)
{
  return(log1p(data*(seq(length(data))*bin_width+range_offset)^2))
}

extinction_coefficient <- function(data,scan_type=NULL,angle=NULL,latest_safe_angle=NULL,latest_safe_measurement=NULL,bin_width,k=1,verbose=FALSE)
{
  #calculating the extinction coefficient from a set of measurements according to Klett inversion
  
  #data preparation
  corrected_data <- range_correction(background_subtraction(data=data,unsupervised=TRUE,bin_width=bin_width,zeroing=FALSE,verbose=verbose),bin_width=bin_width)
  
  #Klett and slope methods were compared during testing - slope method will be removed from the final product
  if (TRUE) #TRUE for Klett inversion, FALSE for slope method
  {
    #in case of atmospheric anomalies, the reference point used in Klett inversion is moved accordingly
    #anomalies <- anomaly_detection(data=expm1(corrected_data),bin_width=bin_width,angle=angle,verbose=verbose)
    anomalies <- anomaly_detection(data=expm1(range_correction(background_subtraction(data=data,unsupervised=TRUE,bin_width=bin_width,zeroing=TRUE,verbose=verbose),bin_width=bin_width)),bin_width=bin_width,angle=angle,verbose=verbose)
    if(sum(abs(anomalies))==0)
      anomalies[length(anomalies)] <- 1
    cutoff <- max(min(seq(length(data))[anomalies==1][1],upper_atmosphere_cutoff(data,bin_width=bin_width)),floor(4500/bin_width))
    
    #reference point as a message
    if (verbose)
    {
      message("Current angle: ",angle)
      message("Optimal reference range after accounting for atmospheric anomalies and a minimum of 4.5 km: ",cutoff*bin_width," m.")
      message("Resulting optimal reference height: ",cutoff*bin_width*sinpi(angle/180)," m.")
    }
    #initialization of reference point
    if(angle == 0)
    {
      latest_safe_angle <- NULL
      latest_safe_measurement <- NULL
    }
    
    if (!(is.null(latest_safe_angle) || is.null(latest_safe_measurement)))
    {
      #if there is a previous "safe measurement, we check whether the current measurement has its reference point below a reasonable height
      if (cutoff*bin_width*sinpi(angle/180) > 3000)
      {
        if (verbose)
        {
          message("Reference point height is adequate.")
          message("Updating safe angle and measurement.")
        }
        #here we update the "safe" measurements
        latest_safe_angle <- angle
        latest_safe_measurement <- data
        #message("Latest safe angle: ",latest_safe_angle)
      } else {
        #if the reference point is too low, we refer to the last "safe" measurement at that height
        if (verbose)
        {
          message("Reference point height too low.")
          message("Latest safe angle: ",latest_safe_angle)
          message("Retrieving latest safe extinction coefficient at a height of ",cutoff*bin_width*sinpi(angle/180)," m.")
          #message("Bin number: ",ceiling(cutoff*sinpi(angle/180)/sinpi(latest_safe_angle/180)))
          message("Closest available extinction coefficient at a height of ", ceiling(cutoff*sinpi(angle/180)/sinpi(latest_safe_angle/180)) * bin_width* sinpi(latest_safe_angle/180)," m.")
          #message("Proxy for local extinction coefficient: ",latest_safe_measurement[ceiling(cutoff*sinpi(angle/180)/sinpi(latest_safe_angle/180))])
        }
      }
      
      #given a reference point, the starting extinction coefficient is calculated via the splope method
      indices <- (cutoff:(cutoff+199))[corrected_data[cutoff:(cutoff+199)]!=0]
      if (length(indices)==0)
        indices <- (cutoff:(cutoff+199))[corrected_data[cutoff:(cutoff+199)]!=0]
      #initialization in case of extremely low signal
      if (length(indices)>0)
      {
        signal <- corrected_data[indices]
        if (length(indices)==1)
        {
          signal <- c(signal,signal*1.001)
          indices <- c(indices,indices+1)
        }
        
        #discarding negative extinction coefficient estimates resulting from noisy data
        #slope method is notoriously unstable at low signal intensities
        if (verbose)
        {
          coefficient_estimate <- (signal[1:(length(signal)-1)]-signal[2:length(signal)])/(indices[2:length(indices)]-indices[1:(length(indices)-1)])/bin_width/2
          if (sum(coefficient_estimate>0)>0)
          {
            coefficient_estimate <- mean(coefficient_estimate[coefficient_estimate>0])
          } else {
            coefficient_estimate <- abs(mean(coefficient_estimate))
          }
          #message("Coefficient estimate: ",coefficient_estimate)
        }
        
        #bimodal coefficient estimate - retrieving the extinction coefficient from a previous measurement if the reference point is below a reasonable height
        #will be used as the default at the final product, after more testing
        if (FALSE || cutoff*bin_width*sinpi(angle/180) > 3000) #FALSE for bimodal calculation
        {
          coefficient_estimate <- (signal[1:(length(signal)-1)]-signal[2:length(signal)])/(indices[2:length(indices)]-indices[1:(length(indices)-1)])/bin_width/2
          if (sum(coefficient_estimate>0)>0)
          {
            coefficient_estimate <- mean(coefficient_estimate[coefficient_estimate>0])
          } else {
            coefficient_estimate <- abs(mean(coefficient_estimate))
          }
        } else {
          if (verbose)
            message("Using proxy for the extinction coefficient instead")
          coefficient_estimate <- latest_safe_measurement[ceiling(cutoff*sinpi(angle/180)/sinpi(latest_safe_angle/180))]
        }
        
        #Klett inversion
        extinction <- exp((corrected_data[min(seq(length(anomalies))[anomalies==1][1],(seq(50)[data[1:50]==max(data[1:50])])[1]):length(corrected_data)]-mean(signal))/k)
        extinction[extinction<0] <- 0
        integrals <- c()
        for (i in 1:length(extinction))
          integrals <- c(integrals,sum(extinction[i:cutoff])*bin_width)
        
        #alternate Klett inversion, correcting for low-signal artifacts concerning the reference point
        temp_extinction <- c(rep(0,min(seq(length(anomalies))[anomalies==1][1],(seq(50)[data[1:50]==max(data[1:50])])[1])-1),extinction/(coefficient_estimate^(-1)+2*integrals/k))
        #coefficient_estimate <- coefficient_estimate * mean(temp_extinction) / max(temp_extinction)
        #coefficient_estimate <- coefficient_estimate /100
        if (log10(max(temp_extinction/mean(temp_extinction))) >= 1.5)
          coefficient_estimate <- coefficient_estimate * 10 ^ ceiling(log10(mean(temp_extinction) / max(temp_extinction)))
        #coefficient_estimate <- coefficient_estimate/100
        #message("Force to local mean coefficient: ", mean(temp_extinction) / max(temp_extinction))
        #message("Order of magnitude correction coefficient: ", 10 ^ ceiling(log10(mean(temp_extinction) / max(temp_extinction))))
        rm(temp_extinction)
        
        #updating the lastest "safe" measurement
        if (cutoff*bin_width*sinpi(angle/180) > 3000)
        {
          if (verbose)
            message("Updating lastest safe measurement.")
          latest_safe_measurement <- c(rep(0,min(seq(length(anomalies))[anomalies==1][1],(seq(50)[data[1:50]==max(data[1:50])])[1])-1),extinction/(coefficient_estimate^(-1)+2*integrals/k))
        }
        return(list(c(rep(0,min(seq(length(anomalies))[anomalies==1][1],(seq(50)[data[1:50]==max(data[1:50])])[1])-1),extinction/(coefficient_estimate^(-1)+2*integrals/k)),latest_safe_angle,latest_safe_measurement))
      } else {
        message("Measurement set contains insufficient non-zero measurements for analysis.")
        return(list(rep(0,length(data)),latest_safe_angle,latest_safe_measurement))
      }
    } else {
      #if there are no previous "safe" measurements, calculation proceeds as normal
      indices <- ((cutoff-199):cutoff)[corrected_data[(cutoff-199):cutoff]!=0]
      if (length(indices)==0)
        indices <- ((cutoff-299):cutoff)[corrected_data[(cutoff-299):cutoff]!=0]
      if (length(indices)>0)
      {
        signal <- corrected_data[indices]
        if (length(indices)==1)
        {
          signal <- c(signal,signal*1.001)
          indices <- c(indices,indices+1)
        }
        coefficient_estimate <- (signal[1:(length(signal)-1)]-signal[2:length(signal)])/(indices[2:length(indices)]-indices[1:(length(indices)-1)])/bin_width/2
        if (sum(coefficient_estimate>0)>0)
        {
          coefficient_estimate <- mean(coefficient_estimate[coefficient_estimate>0])
        } else {
          coefficient_estimate <- abs(mean(coefficient_estimate))
        }
        extinction <- exp((corrected_data[min(seq(length(anomalies))[anomalies==1][1],(seq(50)[data[1:50]==max(data[1:50])])[1]):length(corrected_data)]-mean(signal))/k)
        extinction[extinction<0] <- 0
        integrals <- c()
        for (i in 1:length(extinction))
          integrals <- c(integrals,sum(extinction[i:cutoff])*bin_width)
        
        #alternate Klett inversion, correcting for low-signal artifacts concerning the reference point
        temp_extinction <- c(rep(0,min(seq(length(anomalies))[anomalies==1][1],(seq(50)[data[1:50]==max(data[1:50])])[1])-1),extinction/(coefficient_estimate^(-1)+2*integrals/k))
        #coefficient_estimate <- coefficient_estimate * mean(temp_extinction) / max(temp_extinction)
        #coefficient_estimate <- coefficient_estimate /100
        if (log10(max(temp_extinction/mean(temp_extinction))) >= 1.5)
          coefficient_estimate <- coefficient_estimate * 10 ^ ceiling(log10(mean(temp_extinction) / max(temp_extinction)))
        #coefficient_estimate <- coefficient_estimate/100
        #message("Force to local mean coefficient: ", mean(temp_extinction) / max(temp_extinction))
        #message("Order of magnitude correction coefficient: ", 10 ^ ceiling(log10(mean(temp_extinction) / max(temp_extinction))))
        rm(temp_extinction)
        
        return(list(c(rep(0,min(seq(length(anomalies))[anomalies==1][1],(seq(50)[data[1:50]==max(data[1:50])])[1])-1),extinction/(coefficient_estimate^(-1)+2*integrals/k)),latest_safe_angle,latest_safe_measurement))
      } else {
        message("Measurement set contains insufficient non-zero measurements for analysis.")
        return(list(rep(0,length(data)),latest_safe_angle,latest_safe_measurement))           
      }
    }
  } else {
    #Slope method for extinction coefficient calculation [will be removed from the final product]
    
    #extinction <- as.numeric(na.omit(filter(corrected_data[seq(length(corrected_data)) >= (seq(50)[data[1:50]==max(data[1:50])])[1] & seq(length(corrected_data)) <= min(seq(length(corrected_data))[corrected_data==0])],rep(1,5)/5)))
    #ugly stuff - this was a first attempt
    #extinction_2 <- as.numeric(na.omit(filter(corrected_data,rep(1,5)/5)))[seq(length(as.numeric(na.omit(filter(corrected_data,rep(1,5)/5))))) >= (seq(50)[data[1:50]==max(data[1:50])])[1] & seq(length(as.numeric(na.omit(filter(corrected_data,rep(1,5)/5))))) <= min(seq(length(as.numeric(na.omit(filter(corrected_data,rep(1,5)/5)))))[as.numeric(na.omit(filter(corrected_data,rep(1,5)/5)))==0])]
    #extinction <- c(rep(0,((seq(50)[data[1:50]==max(data[1:50])])[1]-1)[1]),(extinction[2:length(extinction)]-extinction[1:(length(extinction)-1)])/(-2*bin_width*extinction[1:(length(extinction)-1)]),rep(0,length(corrected_data)-(seq(50)[data[1:50]==max(data[1:50])])[1]-length(extinction)+2))
    #return(list(c(rep(0,min(seq(length(anomalies))[anomalies==1][1],(seq(50)[data[1:50]==max(data[1:50])])[1])-1),extinction/(coefficient_estimate^(-1)+2*integrals/k)),latest_safe_angle,latest_safe_measurement))
    #extinction <- extinction
    #extinction <- (corrected_data[((seq(50)[data[1:50]==max(data[1:50])])[1]+1):length(corrected_data)] - corrected_data[(seq(50)[data[1:50]==max(data[1:50])])[1]:(length(corrected_data)-1)])/(2*corrected_data[(seq(50)[data[1:50]==max(data[1:50])])[1]:(length(corrected_data)-1)]*bin_width)
    #barplot(extinction[1:(length(extinction)-1)])
    
    #######second try - every non-zero value
    extinction <- corrected_data[seq(length(corrected_data)) >= (seq(50)[data[1:50]==max(data[1:50])])[1]]
    extinction[length(extinction)] <- 0
    extinction[extinction>0] <- (extinction[seq(length(extinction))[extinction>0]+1]-extinction[extinction>0]) / (-2 * extinction[extinction>0] * bin_width)
    extinction[extinction<0] <- 0
    extinction[extinction > 0.99 * max(extinction)] <- 0
    
    #sum(cumsum(extinction*bin_width)<=3)*bin_width
    #message(length(data))
    #message(((seq(50)[data[1:50]==max(data[1:50])])[1]-1)[1])
    #message(length(list(c(rep(0,(seq(50)[data[1:50]==max(data[1:50])]-1)[1]),extinction),latest_safe_angle,latest_safe_measurement)[[1]]))
    
    return(list(c(rep(0,(seq(50)[data[1:50]==max(data[1:50])]-1)[1]),extinction),latest_safe_angle,latest_safe_measurement))
    #list(c(rep(0,(seq(50)[data[1:50]==max(data[1:50])]-1)[1]),extinction,latest_safe_angle,latest_safe_measurement))[[1]]
  }
}

scanning_profile_extinction <- function(scanning_directory,measurements_of_interest,is_scan,scan_type=NULL,k=1,verbose=FALSE,output_file=FALSE)
{
  #will move this in the arguments, to enable toggling the cloud detection algorithm
  detect_clouds <- TRUE
  
  #cautionary message in case of unprecedented signal variation close to the lidar
  snr_limit <- 0.05
  
  #data acquisition of the first measurement to serve as an initialization template for the data matrix
  file_list <- list.files(scanning_directory)[list.files(scanning_directory)!="Output"]
  temp <- reading_and_preparation(filename=paste(scanning_directory,file_list[1],sep="/"),significant_data_sets=max(measurements_of_interest,2))[,measurements_of_interest]
  
  #initialization of the list of acceptable measurements, giving special consideration to using a single photon counting channel
  if (sum(measurements_of_interest==2) == 1 || sum(measurements_of_interest==4) == 1)
  {
    measurement_check <- rep(sum(reading_and_preparation(filename=paste(scanning_directory,file_list[1],sep="/"),significant_data_sets=2)[,2]!=0)/length(temp) >= 0.05,length(measurements_of_interest))
  } else {
    measurement_check <- rep(TRUE,length(measurements_of_interest))
  }
  
  #data matrix initialization using the first measurement
  if (length(measurements_of_interest)==1)
  {
    data <- matrix(ncol=1,nrow=length(temp))
    data[,1] <- temp
  } else {
    data <- matrix(ncol=length(measurements_of_interest),nrow=dim(temp)[1])
    data[,1:length(measurements_of_interest)] <- temp
  }
  
  #adding the remainder of measurements in the matrix
  bin_width <- as.numeric(read.table(paste(scanning_directory,file_list[1],sep="/"),skip=4,nrows=1,fill=TRUE)[7])
  wavelength <- as.integer(strsplit(as.character(read.table(paste(scanning_directory,file_list[1],sep="/"),skip=4,nrows=1,fill=TRUE)[8]),split=".p")[[1]][1])
  for (i in 2:length(file_list))
  {
    data <- cbind(data,reading_and_preparation(filename=paste(scanning_directory,file_list[i],sep="/"),significant_data_sets=max(measurements_of_interest,2))[,measurements_of_interest])
  
    #initialization of the list of acceptable measurements, giving special consideration to using a single photon counting channel
    if (sum(measurements_of_interest==2) == 1 || sum(measurements_of_interest==4) == 1)
    {
      measurement_check <- c(measurement_check,rep(sum(reading_and_preparation(filename=paste(scanning_directory,file_list[i],sep="/"),significant_data_sets=2)[,2]!=0)/dim(data)[1] >= 0.05,length(measurements_of_interest)))
    } else {
      measurement_check <- c(measurement_check,rep(TRUE,length(measurements_of_interest)))
    }    
  }
  
  #checking for unacceptable measurements
  if (sum(!measurement_check)>0)
  {
    if (sum(measurement_check)>0)
    {
      for (i in 1:(sum(!measurement_check)/length(measurements_of_interest)))
        message("WARNING! Measurement set ",file_list[!measurement_check[seq(length(measurement_check))%%length(measurements_of_interest)==1]][i]," has less than 5% non-zero measurements. Discarded from analysis.")
      data <- data[,measurement_check]
    } else {
      write.table(c("An acceptable measurement set had non-estimable aerosol-free signal. Unable to proceed to analysis.","Heavy visibility obstruction. Visibility less than 100 m, non-estimable.")[as.integer(max(colSums(data[ceiling(200/bin_width):dim(data)[1],]!=0)/dim(data)[1])<=0.01)+1],file="Visibility_unavailable.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
      stop(c("An acceptable measurement set had non-estimable aerosol-free signal. Unable to proceed to analysis.","Heavy visibility obstruction. Visibility less than 100 m, non-estimable.")[as.integer(max(colSums(data[ceiling(200/bin_width):dim(data)[1],]!=0)/dim(data)[1])<=0.01)+1])
    }
  }
  message("Measurements of interest length: ", length(measurements_of_interest))
  if (length(measurements_of_interest) > 1)
    measurement_check <- measurement_check[seq(length(measurement_check))%%length(measurements_of_interest)==1]
  #cat(measurement_check)
  
  if (!is_scan)
  {
    message("Processing data as separate measurements.")
  } else {
    for (i in 1:length(measurements_of_interest))
    {
      #notifying in case of large signal variation close to the lidar
      indices_under_scrutiny <- seq(dim(data)[2])[(seq(dim(data)[2])-i)%%length(measurements_of_interest)==0]
      if (sd(data[1,indices_under_scrutiny]*bin_width^2)/mean(data[1,indices_under_scrutiny]*bin_width^2)>snr_limit)
      {
        message("Measurements in channel ",measurements_of_interest[i]," at a distance of ",bin_width," m exhibit SNR less than ",ceiling(1/snr_limit),".")
        message("Current SNR of channel ",measurements_of_interest[i],": ",floor(mean(data[1,indices_under_scrutiny]*bin_width^2)/sd(data[1,indices_under_scrutiny]*bin_width^2)),". Caution and possible instrument calibrartion recommended.")
      }
    }
  }
  
  #preparing headers according to measurement metadata
  headers <- c()
  for (i in 1:length(file_list[measurement_check]))
    headers <- c(headers,paste(rep(paste("Elevation",-as.numeric(read.table(paste(scanning_directory,file_list[measurement_check][i],sep="/"),skip=1,nrows=1,fill=TRUE)[9]),"Azimuth",(as.numeric(read.table(paste(scanning_directory,file_list[measurement_check][i],sep="/"),skip=1,nrows=1,fill=TRUE)[10])+as.numeric(read.table(paste(scanning_directory,file_list[i],sep="/"),skip=2,nrows=1,fill=TRUE))[7])%%360,sep="_"),length(measurements_of_interest)),measurements_of_interest,sep="_"))
  
  #implementation of total signal calculation if only the appropriate channels are concerned (TRUE to perform)
  if (TRUE && (sum(measurements_of_interest!=c(1,3))==0 || sum(measurements_of_interest!=c(2,4))==0))
  {
    headers <- headers[seq(dim(data)[2])%%2==1]
    data <- data[,seq(dim(data)[2])%%2==1] + 0.06 * data[,seq(dim(data)[2])%%2==0]
    measurements_of_interest <- measurements_of_interest[1]
  }
  
  #if the measurements are handles separately and not as a scan, the filenames are used as headers
  if (is_scan)
  {
    colnames(data) <- headers
  } else {
    colnames(data) <- file_list[measurement_check]
  }
  
  #for (i in 1:length(file_list[measurement_check]))
  #  if (as.numeric(sapply(strsplit(headers,split="_"),'[',2))[i] >= 45)
  #    cloud_detection(filename=colnames(data)[i],data=data[,i],angle=as.numeric(sapply(strsplit(headers,split="_"),'[',2))[i],bin_width,verbose,output_file)
  
  #cloud detection module
  if (detect_clouds)
  {
    if (verbose)
    {
      message("Performing cloud detection.")
    } else {      
      progress <- txtProgressBar(max=length(file_list[measurement_check]),char="=",style=3)
    }
    for (i in 1:length(file_list[measurement_check]))
    {
      if (!verbose)
        setTxtProgressBar(progress,i)
      if (as.numeric(sapply(strsplit(headers,split="_"),'[',2))[i] >= 5)
        cloud_detection(filename=colnames(data)[i],data=data[,i],angle=as.numeric(sapply(strsplit(headers,split="_"),'[',2))[i],bin_width,verbose,output_file)
    }
    if (verbose)
    {
      message("Calculating extinction.")
    } else {
      close(progress)      
    }
  }
  
  
  safe_measurement <- NULL
  safe_angle <- NULL
  
  if (is_scan)
  {
    if (scan_type=="elevation")
    {
      #initialization of the highest "safe" measuremets
      if (verbose)
        message("Scan at ",sapply(strsplit(headers,split="_"),'[',2)[sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix[1]]," degrees can safely used as reference.")
      #message("Initialization.")
      #safe_measurement <- data[,sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix[1]]
      safe_angle <- sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=FALSE)[1]
      safe_measurement <- extinction_coefficient(data=data[,sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix[1]],scan_type=NULL,angle=safe_angle,latest_safe_angle=NULL,latest_safe_measurement=NULL,bin_width=bin_width,k=k,verbose=FALSE)[[1]]
      
      #message("Safe angle: ",safe_angle)
      #message(length(safe_measurement))
      
      #message("Reordering in progress")
      #message(sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix)
      
      #reordering by descending elevation angle
      data <- data[,sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix]
      measurement_check <- measurement_check[sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix]
    } else {
      #in horizontal scans, we do check whether the set includes at least one high enough (this part could be removed from the final version)
      if (!is_scan || max(as.numeric(sapply(strsplit(headers,split="_"),'[',2))) <= 85)
      {
        if (verbose)
          message("No reference scan available, processing low scans individually.")
        safe_angle <- NULL
        safe_measurement <- NULL
      } else {
        if (verbose)
        {
          message("Scan at ",sapply(strsplit(headers,split="_"),'[',2)[sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix[1]]," degrees can safely used as reference.")
          #message("Reference scan measurement at ",sapply(strsplit(headers,split="_"),'[',2)[sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix[1]]," degrees.")
        }
        
        #this part could be removed, too, since azimuth scans do not usually include a vertical one
        safe_angle <- sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=FALSE)[1]
        safe_measurement <- data[,sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix[1]]
        #message("Safe angle: ",safe_angle)
        #message(length(safe_measurement))
      }
    }
  }
  
  temp_safe_measurement <- NULL
  temp_safe_angle <- NULL
  
  if (!verbose)
    progress <- txtProgressBar(max=dim(data)[2],char="=",style=3)
  for (i in 1:dim(data)[2])
  {
    if (verbose && (i%%length(measurements_of_interest)==1 || length(measurements_of_interest)==1))
      message(ceiling(i/length(measurements_of_interest)),"/",sum(measurement_check)," : Processing data file ",file_list[measurement_check][sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix][ceiling(i/length(measurements_of_interest))])
    
    #initialization of "safe" measurement
    if (i==1 && is_scan)
    {
      if (verbose)
        message("Initialization")
      temp_safe_measurement <- safe_measurement
      temp_safe_angle <- safe_angle
    }
    if (verbose && is_scan && !is.null(temp_safe_angle))
      message("Lowest safe angle: ",temp_safe_angle)
    
    #updating the "safe" measurement
    if ((as.integer(strsplit(headers[sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix][i],split="_")[[1]][2])) > 0)
    {
      temp_data <- extinction_coefficient(data=data[,i],scan_type=scan_type,angle=as.integer(strsplit(headers[sort(as.numeric(sapply(strsplit(headers,split="_"),'[',2)),decreasing=TRUE,index.return=TRUE)$ix][i],split="_")[[1]][2]),latest_safe_angle=temp_safe_angle,latest_safe_measurement=temp_safe_measurement,bin_width=bin_width,k=k,verbose=verbose)
      
      temp_safe_angle <- temp_data[[2]]
      temp_safe_measurement <- temp_data[[3]]
      temp_data <- temp_data[[1]] #needs tidying up
      data[,i] <- temp_data
      #message("Safe angle to be used next: ",temp_safe_angle)
    }
    
    if (!verbose)
      setTxtProgressBar(progress,i)
  }
  rm(safe_angle,temp_safe_angle,safe_measurement,temp_safe_measurement)
  
  #measurement_check[measurement_check][colSums(data)==0] <- FALSE
  data <- data[,colSums(data)!=0]
  if (!verbose)
    close(progress)
  rownames(data) <- as.character(seq(dim(data)[1])*bin_width)
  if(output_file)
    write.table(data,file=paste("Radial_extinction_coefficients",wavelength,"nm.txt",sep="_"),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
  #data <- data[,colSums(data*bin_width)>=3]
  return(data)
}

visibility_range <- function(extinction,bin_width,model=NULL,wavelength,incoming=FALSE,incoming_range=NULL,verbose=FALSE)
{
  #this method calculates the meteorological visibility according to the Koschmieder formula and ICAO argument for 2% threshold
  
  if (incoming && is.null(incoming_range))
    stop("Please provide the distance of the incoming object in metres.")
  if (verbose)
  {
    if (!is.null(model))
    {
      message("Selected model: ",as.character(model))
      message("Converting visibility from ",wavelength," nm to the visible spectrum.")
    } else {
      message("No model selected. Visibility will be calculated at ",wavelength," nm.")
    }
  }
  
  #if the system wavelength has an associated empirical conversion-to-human-range model, it is implemented, otherwise the Angstrom Exponent method is used
  if (wavelength==355)
  {
    if (!is.null(model))
      extinction <- angstrom_exponent_extinction_coefficient_conversion(extinction,wavelength)
  } else if (wavelength==1064) {
    if (!is.null(model))
    {
      if (model=="urban-rural")
        extinction <- (extinction/0.503)^(1/1.08)
      if (model=="maritime")
        extinction <- (extinction/0.983)^(1/1.16)
      if (model=="angstrom_exponent")
        extinction <- angstrom_exponent_extinction_coefficient_conversion(extinction,wavelength)
    }
  } else if (wavelength==1550 | wavelength==1570) {
    if (!is.null(model))
    {
      if (model=="urban-rural")
        extinction <- (extinction/0.314)^(1/1.11)
      if (model=="maritime")
        extinction <- (extinction/0.996)^(1/1.20)
      if (model=="angstrom_exponent")
        extinction <- angstrom_exponent_extinction_coefficient_conversion(extinction,wavelength)
    }
  }
  
  #implementation of the Koschmieder formula, returning visibility range im metres
  visibility <- c()
  if (incoming)
  {
    for (i in floor(incoming_range/bin_width):1)
      visibility <- c(visibility,sum(extinction[floor(incoming_range/bin_width):i]*bin_width))
    if (visibility[length(visibility)]>=3)
    {
      return(sum(visibility<=3)*bin_width)
    } else {
      return(c(incoming_range,visibility[length(visibility)]))
    }
  } else {
    for (i in 1:length(extinction))
    {
      visibility <- c(visibility,sum(extinction[1:i]*bin_width))
      if (visibility[i]>3)
        break
    }
    return(sum(visibility<=3)*bin_width)
  }
}

progressive_slant_range <- function(cartesian_extinction,bin_width,model=NULL,wavelength,incoming=FALSE,incoming_distance,incoming_height,verbose=FALSE)
{
  #the aim is to find a reasonable slant range by starting from the ground direcly below and increasing the angle until the MOR becomes less than the distance to the ground at that angle
  #in all cases, in the interests of speed, the angle is calculated within 45 degrees, then within 15 degrees, then within 5 degrees and lastly within 1 degree
  
  #reversed reflects whether the search has exceeded 45 degrees in angle
  reversed <- FALSE
  
  #visible_airport refers to whether an incoming craft has a slant range long enough to see the airport. If so, the slant range is set to the distance to the airport.
  visible_airport <- FALSE
  
  #bahroken refers to the endge cases where the slant range reaches the geometric limits of available data 
  bahroken <- FALSE
  
  if (incoming)
  {
    #checking whether the airport can be seen in the case of an incoming aircraft and setting the slant range accordingly
    #message("Incoming, checking whether airport is visible")
    
    if (incoming_height==incoming_distance)
    {
      #message("Airport at 45 degrees")
      angle <- 45
      progressive_slant_visibility_temp <- visibility_range(extinction=diag(cartesian_extinction[1:ceiling(incoming_height/bin_width),1:ceiling(incoming_distance/bin_width)]),bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)
    } else {
      if (incoming_distance > incoming_height)
      {
        #message("Airport at more than 45 degrees")
        angle <- atan(incoming_height/incoming_distance)/pi*180
        minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),1:ceiling(incoming_distance/bin_width)])
        offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
        
        angled_extinction_profile <- c()
        for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
        {
          angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
          #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
        }
      } else {
        temp <- incoming_height
        incoming_height <- incoming_distance
        incoming_distance <- temp
        reversed <- TRUE
        angle <- atan(incoming_height/incoming_distance)/pi*180
        minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),1:ceiling(incoming_distance/bin_width)])
        offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
        
        angled_extinction_profile <- c()
        for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
        {
          angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
          #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
        }
      }
      
      progressive_slant_visibility_temp <- visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)
      if (reversed)
      {
        temp <- incoming_height
        incoming_height <- incoming_distance
        incoming_distance <- temp
      }
    }
    
    #given that the airport is beyond the slant range, the angle is progressively narrowed down
    if (is.na(progressive_slant_visibility_temp[2]))
    {
      #message("Can't see the airport")
      #the airport is not visible, and the edge of the data is at less than 45 degrees
      if (incoming_distance <= incoming_height) 
      {
        #message("Airport invisible and at less than 45 degrees from vertical")
        #message("Therefore, searching at less than 45 degrees")
        
        #less than 45 degrees
        for (i in 1:2)
        {
          angle <- 15*i
          if (angle >= atan(incoming_distance/incoming_height)/pi*180)
          {
            true_angle <- angle
            break
          }
          minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
          
          angled_extinction_profile <- c()
          offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
          for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
          {
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
          }
          #message(angle)
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
            break
          angle <- 45
        }
        true_angle <- angle-15
        
        #more than 0, 15 or 30 and less than 45
        for (i in 1:2)
        {
          angle <- true_angle + 5*i
          if (angle >= atan(incoming_distance/incoming_height)/pi*180)
          {
            true_angle <- angle
            break
          }
          minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
          
          angled_extinction_profile <- c()
          offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
          for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
          {
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
          }
          #message(angle)
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
            break
          angle <- true_angle + 15
        }
        
        true_angle <- angle-5
        
        #last 5
        for (i in 1:4)
        {
          angle <- true_angle + i
          if (angle >= atan(incoming_distance/incoming_height)/pi*180)
          {
            true_angle <- angle
            break
          }
          minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
          
          angled_extinction_profile <- c()
          offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
          for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
          {
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
          }
          #message(angle)
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
            break
          angle <- true_angle + 5
        }
        true_angle <- angle-1 #this is the final angle
        
      } else {        
        #the airport is not visible, and the edge of the data is at above 45 degrees
        
        #message("Airport invisible and at more than 45 degrees from vertical")
        #airport is not visible, and the edge is above 45 degrees
        #therefore, searching below and above the diagonal
        angle <- 45
        if (is.na(visibility_range(extinction=diag(cartesian_extinction[1:ceiling(incoming_height/bin_width),seq(sort(c(ceiling(incoming_distance/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * ceiling(incoming_height/bin_width)))[1] +1,sort(c(ceiling(incoming_distance/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * ceiling(incoming_height/bin_width)))[2])]),bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
        {
          #message("Can't see the ground at 45 degrees")
          #less than diagonal
          for (i in 1:2)
          {
            angle <- 15*i
            minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
            
            angled_extinction_profile <- c()
            offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
            for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
            {
              angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
              #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
            }
            #message(angle)
            if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
              break
            #angle <- 15
          }
          true_angle <- angle-15
          
          #more than 0, 15 or 30 and less than 45
          for (i in 1:2)
          {
            angle <- true_angle + 5*i
            minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
            
            angled_extinction_profile <- c()
            offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
            for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
            {
              angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
              #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
            }
            #message(angle)
            if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
              break
            angle <- true_angle + 15
          }
          
          true_angle <- angle-5
          
          #last 5
          for (i in 1:4)
          {
            angle <- true_angle + i
            minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
            
            angled_extinction_profile <- c()
            offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
            for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
            {
              angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
              #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
            }
            #message(angle)
            if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
              break
            angle <- true_angle + 5
          }
          true_angle <- angle-1 #this is the final angle
        } else {
          #airport is beyond the slant range, true angle is over 45 degrees, the aircraft may be incoming or outcoming
          #message("Can see the ground at 45 degrees")
          
          for (i in 2:1)
          {
            angle <- 15*i
            if (c(incoming_distance,dim(cartesian_extinction)[2]*bin_width-incoming_distance)[as.integer(!incoming)+1] <= incoming_height/tan((angle)*pi/180))
              break
            #minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]])
            #minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),c(1,ceiling(incoming_height/bin_width))[as.integer(!incoming)+1]:c(ceiling(incoming_distance/bin_width),dim(cartesian_extinction)[2])[as.integer(!incoming)+1]])
            ##########################################################
            #range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]
            minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))])
            angled_extinction_profile <- c()
            offsets <- ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))
            for (j in range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[2])
            {
              angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j-range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]+1)])
              #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
            }
            #message(angle)
            #message(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2])
            if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2]))
              break
            angle <- 0
          }
          true_angle <- angle
          
          #less than 0, 15 or 30 and less than 45
          for (i in 2:1)
          {
            angle <- true_angle + 5*i
            if (c(incoming_distance,dim(cartesian_extinction)[2]*bin_width-incoming_distance)[as.integer(!incoming)+1] <= incoming_height/tan((angle)*pi/180))
              break
            minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))])
            angled_extinction_profile <- c()
            offsets <- ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))
            
            for (j in range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[2])
            {
              angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j-range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]+1)])
              #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
            }
            #message(angle)
            #message(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2])
            if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2]))
              break
            angle <- true_angle
          }
          true_angle <- angle
          
          #last 5
          for (i in 4:1)
          {
            angle <- true_angle + i
            minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))])
            angled_extinction_profile <- c()
            offsets <- ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))
            
            for (j in range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[2])
            {
              angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j-range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]+1)])
              #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
            }
            #message(angle)
            #message(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2])
            if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2]))
              break
            angle <- true_angle
          }
          true_angle <- 90-angle #this is the final angle
          if (c(incoming_distance,dim(cartesian_extinction)[2]*bin_width-incoming_distance)[as.integer(!incoming)+1] <= incoming_height/tan((angle)*pi/180))
          {
            visible_airport <- TRUE
            #break
          }
        }
      }
    } else {
      #message("Can see the airport")
      visible_airport <- TRUE
    }
  } else {
    #message("Outcoming")
    #beyond this point, we deal exclusively with outcoming aircraft, so there is no need to limit the slant range to the distance to the airport
    
    angle <- 45
    if (!c(incoming_distance,dim(cartesian_extinction)[2]*bin_width-incoming_distance)[as.integer(!incoming)+1] <= incoming_height*tan((angle)*pi/180))
    {
      #checking whether the aircarft is near the far edge of the data
      #message("Can't see the end of the map at 45 degrees")
      if (is.na(visibility_range(extinction=diag(cartesian_extinction[1:ceiling(incoming_height/bin_width),seq(sort(c(ceiling(incoming_distance/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * ceiling(incoming_height/bin_width)))[1] +1,sort(c(ceiling(incoming_distance/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * ceiling(incoming_height/bin_width)))[2])]),bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
      {
        #message("Can't see the ground at 45 degrees")
        #less than diagonal
        for (i in 1:2)
        {
          angle <- 15*i
          minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
          angled_extinction_profile <- c()
          offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
          for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
          {
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
          }
          #message(angle)
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
            break
          #angle <- 15
        }
        true_angle <- angle-15
        
        #more than 0, 15 or 30 and less than 45
        for (i in 1:2)
        {
          angle <- true_angle + 5*i
          minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
          angled_extinction_profile <- c()
          offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
          for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
          {
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
            
          }
          #message(angle)
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
            break
          angle <- true_angle + 15
        }
        
        true_angle <- angle-5
        
        #last 5
        for (i in 1:4)
        {
          angle <- true_angle + i
          minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
          angled_extinction_profile <- c()
          offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
          for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
          {
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
          }
          #message(angle)
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
            break
          angle <- true_angle + 5
        }
        true_angle <- angle-1 #this is the final angle
      } else {
        #the ground at 45 degrees can be seen, searching accordingly
        #message("Can see the ground at 45 degrees, searching over them")
        for (i in 2:1)
        {
          angle <- 15*i
          #minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]])
          #minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),c(1,ceiling(incoming_height/bin_width))[as.integer(!incoming)+1]:c(ceiling(incoming_distance/bin_width),dim(cartesian_extinction)[2])[as.integer(!incoming)+1]])
          ##########################################################
          #range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]
          minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))])
          angled_extinction_profile <- c()
          offsets <- ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))
          for (j in range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[2])
          {
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j-range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
          }
          #message(angle)
          #message(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2])
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2]))
            break
          angle <- 0
        }
        true_angle <- angle
        
        #less than 0, 15 or 30 and less than 45
        for (i in 2:1)
        {
          angle <- true_angle + 5*i
          minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))])
          angled_extinction_profile <- c()
          offsets <- ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))
          for (j in range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[2])
          {
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j-range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
            
          }
          #message(angle)
          #message(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2])
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2]))
            break
          angle <- true_angle
        }
        true_angle <- angle
        
        #last 5
        for (i in 4:1)
        {
          angle <- true_angle + i
          minnie_profile <- t(cartesian_extinction[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))])
          angled_extinction_profile <- c()
          offsets <- ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width)))
          for (j in range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[2])
          {
            #minnie_beginning <- 1+sum(offsets<j)
            #if (minnie_beginning > dim(minnie_profile)[1])
            #  minnie_beginning <- dim(minnie_profile)[1]
            #minnie_end <- sum(offsets<j+1)
            #if (minnie_end > dim(minnie_profile)[1])
            #  minnie_end <- dim(minnie_profile)[1]
            #minnie_length <- j-range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]+1
            #if (minnie_length > dim(minnie_profile)[2])
            #  minnie_length <- dim(minnie_profile)[2]
            
            angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j-range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(floor(incoming_height/tan(angle*pi/180)/bin_width))%/%ceiling(floor(incoming_height/tan(angle*pi/180)/bin_width)/ceiling(tan(angle*pi/180)*floor(incoming_height/tan(angle*pi/180)/bin_width))))[1]+1)])
            #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
          }
          #message(angle)
          #message(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2])
          if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/sin(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/sin(angle*pi/180),verbose=FALSE)[2]))
            break
          angle <- true_angle
        }
        true_angle <- 90-angle #this is the final angle
        if (c(incoming_distance,dim(cartesian_extinction)[2]*bin_width-incoming_distance)[as.integer(!incoming)+1] <= incoming_height/tan((angle)*pi/180))
        {
          bahroken <- TRUE
          #break
        }
      }
    } else {
      #this is the edge case of an aircraft near the far edge of the data
      #message("Can see the end of the map at 45 degrees, searching in there")
      #less than diagonal
      for (i in 1:2)
      {
        angle <- 15*i
        if (c(incoming_distance,dim(cartesian_extinction)[2]*bin_width-incoming_distance)[as.integer(!incoming)+1] <= incoming_height*tan((90-angle)*pi/180))
          break
        minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
        angled_extinction_profile <- c()
        offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
        for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
        {
          angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
          #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
        }
        #message(angle)
        if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
          break
        #angle <- 15
      }
      true_angle <- angle-15
      
      #more than 0, 15 or 30 and less than 45
      for (i in 1:2)
      {
        angle <- true_angle + 5*i
        if (c(incoming_distance,dim(cartesian_extinction)[2]*bin_width-incoming_distance)[as.integer(!incoming)+1] <= incoming_height*tan((90-angle)*pi/180))
          break
        minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(!incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
        
        angled_extinction_profile <- c()
        offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
        for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
        {
          angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
          #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
        }
        #message(angle)
        if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
          break
        angle <- true_angle + 15
      }
      
      true_angle <- angle-5
      
      #last 5
      for (i in 1:4)
      {
        angle <- true_angle + i
        #message(angle)
        minnie_profile <- cartesian_extinction[1:ceiling(incoming_height/bin_width),range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(ceiling(incoming_distance/bin_width) + c(-1,1)[as.integer(incoming)+1] * seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2]]
        angled_extinction_profile <- c()
        #change in offset
        offsets <- seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width)))
        for (j in range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[1]:range(seq(ceiling(incoming_height/bin_width))%/%ceiling(ceiling(incoming_height/bin_width)/ceiling(tan(angle*pi/180)*ceiling(incoming_height/bin_width))))[2])
        {
          angled_extinction_profile <- c(angled_extinction_profile,minnie_profile[min(dim(minnie_profile)[1],1+sum(offsets<j)):min(dim(minnie_profile)[1],sum(offsets<j+1)),min(dim(minnie_profile)[2],j+1)])
          #angled_extinction_profile <- c(angled_extinction_profile,rep(colnames(minnie_profile)[j+1],1+diff(range((1+sum(offsets<j)):sum(offsets<j+1)))))
          
        }
        #message(angle)
        if (is.na(visibility_range(extinction=angled_extinction_profile,bin_width=bin_width/cos(angle*pi/180),model,wavelength,incoming=TRUE,incoming_range=incoming_height/cos(angle*pi/180),verbose=FALSE)[2]))
          break
        angle <- true_angle + 5
      }
      
      if (c(incoming_distance,dim(cartesian_extinction)[2]*bin_width-incoming_distance)[as.integer(!incoming)+1] <= incoming_height*tan((90-angle)*pi/180))
      {
        bahroken <- TRUE
        #break
      }
      true_angle <- angle-1 #this is the final angle
    }
  }
  
  if (!visible_airport & !bahroken)
  {
    #message("True angle: ",true_angle)
    return(incoming_height*tan(true_angle*pi/180))
  } else {
    if (visible_airport)
    {
      #the aircraft is approaching the airport and can see it
      return(-1)
    } else {
      #the aircraft can see the far edge of the data - here be dragons
      return(-2)
    }
  }
}

angstrom_exponent_extinction_coefficient_conversion <- function(extinction,lidar_wavelength)
{
  #the Angstrom Exponent methodology is used as a default way to convert results obtained beyond the visible spectrum to results meaningful to human sight
  
  #a table of AERONET historical data is used in order to both facilitate the calculation of the exponents AND to enable the usage of localized data in the interests of improved accuracy
  optical_depth_data <- read.table("/home/SAFETRANS/AERONET_data.txt",header=TRUE,sep="\t")
  #optical_depth_data <- read.table("C:/Users/Hector-Xavier/Desktop/Default directory/AERONET_data.txt",header=TRUE,sep="\t")
  optical_depth_data <- optical_depth_data[,!is.na(colSums(optical_depth_data))]
  
  #this part has been written specifically for the wavelengths proposed in the entiredy of the SAFETRANS program
  #different wavelengths will use different rows
  if (lidar_wavelength==355)
  {
    wavelengths <- c(1,2)
  } else if (lidar_wavelength==1064) {
    wavelengths <- c(6,7)
  } else if (lidar_wavelength==1550 | lidar_wavelength==1570) {
    wavelengths <- c(7,6)
  }
  starting_wavelength <- as.integer(rownames(optical_depth_data)[wavelengths[1]])
  if (lidar_wavelength==1064)
    wavelengths <- wavelengths + 1
  lidar_optical_depth <- as.numeric(optical_depth_data[wavelengths[1],])*(lidar_wavelength/starting_wavelength)^(log(as.numeric(optical_depth_data[wavelengths[2],]/optical_depth_data[wavelengths[1],]))/log(as.integer(rownames(optical_depth_data)[wavelengths[2]])/lidar_wavelength))
  wavelengths <- c(4,5)
  visible_optical_depth <- as.numeric(optical_depth_data[wavelengths[1],])*(550/as.integer(rownames(optical_depth_data)[wavelengths[1]]))^(log(as.numeric(optical_depth_data[wavelengths[2],]/optical_depth_data[wavelengths[1],]))/log(as.integer(rownames(optical_depth_data)[wavelengths[2]])/550))
  coefficient <- mean((550/lidar_wavelength)^(log(visible_optical_depth/lidar_optical_depth)/log(550/lidar_wavelength)))
  return(extinction*coefficient)
}

radial_visibility_profile <- function(extinction_profile,is_scan=TRUE,model=NULL,wavelength,output_file=TRUE,verbose=FALSE,research_material=FALSE)
{
  #this is a wrapper than makes use of the other functions to calculate visibility as a set of radial measurements, such as an azimuth scan
  
  bin_width <- as.numeric(rownames(extinction_profile)[2])-as.numeric(rownames(extinction_profile)[1])
  visibility <- c()
  if (verbose)
  {
    for (i in 1:dim(extinction_profile)[2])
    {
      message(i,"/",dim(extinction_profile)[2]," - Measurement set designation: ",colnames(extinction_profile)[i])
      visibility <- c(visibility,visibility_range(extinction=extinction_profile[,i],bin_width,model=model,wavelength=wavelength,verbose=TRUE)[1])
      message("Outward visibility: ",visibility[i]," m.")
    }
  } else {
    progress <- txtProgressBar(max=dim(extinction_profile)[2],char="=",style=3)
    for (i in 1:dim(extinction_profile)[2])
    {
      visibility <- c(visibility,visibility_range(extinction=extinction_profile[,i],bin_width,model=model,wavelength=wavelength)[1])
      setTxtProgressBar(progress,i)
    }
    close(progress)
  }
  if (is.null(model))
    model <- "no_model"
  
  #in the case of requesting output in the form of files, a table of visibilities is produced
  if (output_file)
    if(!file.exists(paste("Radial_outward_visibility_distance_",model,".txt",sep="")))
      write.table(visibility,file=paste("Radial_outward_visibility_distance_",model,".txt",sep=""),quote=FALSE,sep="\t",row.names=colnames(extinction_profile),col.names="Visibility_in_metres")
  visibility <- matrix(visibility,ncol=1,dimnames=list(colnames(extinction_profile),model))
  if (output_file)
  {
    #in the case of requesting output in the form of files, plots are produced
    #(output_file && is_scan)
    extinction_profile <- extinction_profile[,visibility < dim(extinction_profile)[1]*bin_width]
    visibility <- visibility[visibility < dim(extinction_profile)[1]*bin_width,1]
    visibility <- matrix(visibility,ncol=1,dimnames=list(colnames(extinction_profile),model))
    if(!file.exists(paste(getwd(),"Azimuth_visibility_plots",sep="/")))
      dir.create(paste(getwd(),"Azimuth_visibility_plots",sep="/"))
    
    if (is_scan)
    {
      #creation of a a radial visibility plot
      
      angle <- c()
      channels <- c()
      for (i in 1:length(visibility))
      {
        angle <- c(angle,as.integer(strsplit(rownames(visibility)[i],split="_")[[1]][4]))
        channels <- c(channels,as.integer(strsplit(rownames(visibility)[i],split="_")[[1]][5]))
      }
      for (i in 1:length(levels(as.factor(channels))))
      {
        warn_default <- getOption("warn")
        options(warn=-1)
        #png(file=file.path(paste(getwd(),"Azimuth_visibility_plots",sep="/"),paste("Radial_visibility_",model,"_channel_",levels(as.factor(channels))[i],".png", sep = "")),width=1000,height=1000)
        png(file=file.path(paste(getwd(),"Azimuth_visibility_plots",sep="/"),paste("Radial_visibility_",model,"_channel_",levels(as.factor(channels))[i],".png", sep = "")),width=870,height=870)
        #radial.plot(lengths=visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1],radial.pos=angle[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i]/180*pi,labels=c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW"),label.pos=(seq(16)-1)/8*pi,start=+pi/2,clockwise=TRUE,rp.type="pt",point.symbols=visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1],point.col=c("green4","black","darkorange","red")[rowSums(visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1] <= matrix(rep(c(1000,6000,10000),length(visibility[,1])),ncol=3,byrow=TRUE)) + 1],radial.lim=pretty(range(c(0,visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1]))),show.grid.labels=0,radial.labels=paste(pretty(range(c(0,visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1])))/1000,"km",sep=" "),show.centroid=TRUE,main=paste("Radial visibility at",strsplit(rownames(visibility)[1],split="_")[[1]][2],"degrees in channel",levels(as.factor(channels))[i],sep=" "))
        radial.plot(lengths=visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1],radial.pos=angle[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i]/180*pi,labels=c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW"),label.pos=(seq(16)-1)/8*pi,start=+pi/2,clockwise=TRUE,rp.type="pts",point.symbols=visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1],point.col=c("green4","black","darkorange","red")[rowSums(visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1] <= matrix(rep(c(4500,7000,10000),length(visibility[,1])),ncol=3,byrow=TRUE)) + 1],radial.lim=pretty(range(c(0,visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1]))),show.grid.labels=0,radial.labels=paste(pretty(range(c(0,visibility[(seq(length(channels)/length(levels(as.factor(channels))))-1)*length(levels(as.factor(channels)))+i,1])))/1000,"km",sep=" "),show.centroid=TRUE,main=paste("Radial visibility at",strsplit(rownames(visibility)[1],split="_")[[1]][2],"degrees in channel",levels(as.factor(channels))[i],sep=" "))
        legend("topleft",legend=c("Visibility over 10 km","Visibility between 7km and 10km","Visibility between 4.5km and 7km","Visibility less than 4.5km"),pch=15,col=c("green4","black","darkorange","red"))
        dev.off()
        options(warn = warn_default)
        rm(warn_default)
      }
    }
    
    if(research_material && !file.exists(paste(getwd(),"Azimuth_visibility_plots",paste("Incoming_visibilities",model,sep="_"),sep="/")))
    {
      #creation of an "incoming object's visibility barplot" for each measurement, mainly for research purposes
      
      dir.create(paste(getwd(),"Azimuth_visibility_plots",paste("Incoming_visibilities",model,sep="_"),sep="/"))
      message("Calculating radial visibility ranges of incoming objects. Selected atmospheric model: ",model,".") 
      if (!verbose)
        progress <- txtProgressBar(max=length(visibility),char="=",style=3)
      for (i in 1:length(visibility))
      {
        if (verbose)
        {
          message(i,"/",dim(extinction_profile)[2]," - Calculating visibilities for measurement set: ",colnames(extinction_profile)[i])
        } else {
          setTxtProgressBar(progress,i)
        }
        if (model=="no_model")
          model <- NULL
        incoming_vis <- c()
        for (j in seq(6000,1,-50))
          incoming_vis <- c(incoming_vis,visibility_range(extinction_profile[,i],bin_width,model,wavelength,incoming=TRUE,incoming_range=j*bin_width,verbose=FALSE)[1])
        if (is.null(model))
          model <- "no_model"
        png(file=file.path(paste(paste(getwd(),"Azimuth_visibility_plots",paste("Incoming_visibilities",model,sep="_"),sep="/"),paste("Visibility_",rownames(visibility)[i],".png", sep = ""),sep="/")),width=2000,height=1600)
        par(mar= c(6.5, 6.5, 4, 2) + 0.1)
        par(mgp = c(5,1,0))
        barplot(incoming_vis[length(incoming_vis):1]/1000,log="x",col=c("red","grey")[(seq(length(incoming_vis))*50*bin_width>=visibility[i,1])+1],names.arg=paste(seq(6000,1,-50)[length(incoming_vis):1]*bin_width/1000,"km",sep=" "),horiz=TRUE,las=2,xlab="Visibility [km]",ylab="Distance",main=paste(rownames(visibility)[i], " [model: ",model,"]",sep=""))
        legend("topleft",legend=c("Optical contact with overhead airspace","No optical contact"),pch=15,col=c("red","grey"))
        par(mar= c(5, 4, 4, 2) + 0.1)
        par(mgp = c(3,1,0))
        dev.off()
      }
      if (!verbose)
        close(progress)
    } else {
      if(research_material)
        message("Incoming radial visibility plots for atmospheric model '",model,"' already exist. Calculation & visualization skipped.")
    }
  }
  #return(visibility)
}

cartesian_visibility_profile <- function(extinction_profile,model=NULL,wavelength,incoming=FALSE,incoming_distance=NULL,incoming_height=NULL,output_files=TRUE,verbose=TRUE,research_material=FALSE)
{
  #this is a wrapper than makes use of the other functions to calculate visibility as a set derived from a zenith scan
  
  if (incoming && (is.null(incoming_distance) || is.null(incoming_height)))
    stop("Please designate both the height and the distance of the incoming object.")
  maximum_height <- 15000
  bin_width <- as.numeric(rownames(extinction_profile)[2])-as.numeric(rownames(extinction_profile)[1])
  
  #setting up a 2-D cartesian representation of the airspace along the plane of the zenith scan
  cartesian_profile <- matrix(0,nrow=ceiling(maximum_height/bin_width),ncol=dim(extinction_profile)[1])
  colnames(cartesian_profile) <- paste("Distance",rownames(extinction_profile)[1:dim(cartesian_profile)[2]],"m",sep="_")
  rownames(cartesian_profile) <- paste("Height",rownames(extinction_profile)[1:dim(cartesian_profile)[1]],"m",sep="_")
  if (incoming && (ceiling(incoming_distance/bin_width)>dim(cartesian_profile)[2] || incoming_height>maximum_height))
    stop("Designated coordinates beyond the detection range. Maximum height: ",maximum_height," m, maximum distance: ",dim(cartesian_profile)[2]*bin_width," m.")
  
  #filling up the table with the data according to each measurement's angle
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
  
  #filling up the spaces between measurements by assuming horizontally linear changes of the extinction coefficient values
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
  if (sum(rowSums(cartesian_profile)==0)==1)
  {
    if (seq(dim(cartesian_profile)[1])[rowSums(cartesian_profile)==0][1]==dim(cartesian_profile)[1])
    {
      non_zero <- seq(dim(cartesian_profile)[1])[rowSums(cartesian_profile)==0][1]-1
    } else {
      non_zero <- seq(dim(cartesian_profile)[1])[rowSums(cartesian_profile)==0][1]+1
    }
  }
  for (i in seq(dim(cartesian_profile)[1])[rowSums(cartesian_profile)==0][1]:(seq(dim(cartesian_profile)[1])[rowSums(cartesian_profile)==0][sum(rowSums(cartesian_profile)==0)]))
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
  
  #here, the visibility ranges of an incoming aircraft are calculated (horizontal, , vertical and slant)
  if (!is.null(incoming_height) && !is.null(incoming_distance))
  {
    if(incoming)
    {
      horizontal_visibility <- min(visibility_range(extinction=cartesian_profile[ceiling(incoming_height/bin_width),1:ceiling(incoming_distance/bin_width)],bin_width,model,wavelength,incoming,incoming_distance,verbose)[1],incoming_distance)
    } else {
      horizontal_visibility <- min(visibility_range(extinction=c(cartesian_profile[ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width):dim(cartesian_profile)[2]],rep(cartesian_profile[ceiling(incoming_height/bin_width),dim(cartesian_profile)[2]],ceiling(20000/bin_width))),bin_width,model,wavelength,incoming,incoming_distance,verbose)[1],dim(cartesian_profile)[2]*bin_width-incoming_distance)
    }
    vertical_visibility <- visibility_range(extinction=cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)],bin_width,model,wavelength,incoming=TRUE,incoming_height,verbose=FALSE)
    #if(is.null(vertical_visibility[2]))
    if (vertical_visibility[1] < incoming_height)    
    {
      median_slant_visibility <- "No optical contact between object and ground. Median slant visibility unavailable."
      minimum_slant_visibility <- "No optical contact between object and ground. Minimum slant visibility unavailable."
      homogeneous_slant_visibility <- "No optical contact between object and ground. Homogeneous slant visibility unavailable."
      progressive_slant_visibility <- "No optical contact between object and ground. Progressive slant visibility unavailable."
    } else {
      #median slant visibility is intended as a quality control while the development is active (extends the slant range using the median value to get a low slant estimate)
      median_pseudo_visibility <- visibility_range(extinction=c(rep(median(cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),maximum_height/bin_width),cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),bin_width,model,wavelength,incoming=TRUE,incoming_height+maximum_height,verbose=FALSE)[1]
      if (median_pseudo_visibility==incoming_height+maximum_height)
        median_pseudo_visibility <- visibility_range(extinction=c(rep(median(cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),dim(cartesian_profile)[2]),cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),bin_width,model,wavelength,incoming=TRUE,dim(cartesian_profile)[2]*bin_width+incoming_height,verbose=FALSE)[1]
      median_slant_visibility <- ceiling(sqrt(median_pseudo_visibility^2 - incoming_height^2))
      
      #minimum slant visibility is intended as a quality control while the development is active (extends the slant range using the minimum value to get a high slant estimate)
      minimum_pseudo_visibility <- visibility_range(extinction=c(rep(max(cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),maximum_height/bin_width),cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),bin_width,model,wavelength,incoming=TRUE,incoming_height+maximum_height,verbose=FALSE)[1]
      if (minimum_pseudo_visibility==incoming_height+maximum_height)
        minimum_pseudo_visibility <- visibility_range(extinction=c(rep(max(cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),dim(cartesian_profile)[2]),cartesian_profile[1:ceiling(incoming_height/bin_width),ceiling(incoming_distance/bin_width)]),bin_width,model,wavelength,incoming=TRUE,dim(cartesian_profile)[2]*bin_width+incoming_height,verbose=FALSE)[1]
      minimum_slant_visibility <- ceiling(sqrt(minimum_pseudo_visibility^2 - incoming_height^2))
      
      #actual proposed slant range algorithm execution
      progressive_slant_visibility <- progressive_slant_range(cartesian_extinction=cartesian_profile,bin_width,model,wavelength,incoming,incoming_distance,incoming_height,verbose)
      if (progressive_slant_visibility < 0)
        progressive_slant_visibility <- c(incoming_distance,dim(cartesian_profile)[2]*bin_width-incoming_distance)[abs(progressive_slant_visibility)]
      
      #this is the ICAO algorithm for slant visibility, assuming full horizontal homogeinity of the extinction coefficients
      homogeneous_slant_visibility <- floor(incoming_height * tan(acos(vertical_visibility[2]/3)))
      
      if (incoming)
      {
        median_slant_visibility <- min(median_slant_visibility,incoming_distance)
        minimum_slant_visibility <- min(minimum_slant_visibility,incoming_distance)
        homogeneous_slant_visibility <- min(homogeneous_slant_visibility,incoming_distance)
      } else {
        median_slant_visibility <- min(median_slant_visibility,dim(cartesian_profile)[2]*bin_width-incoming_distance)
        minimum_slant_visibility <- min(minimum_slant_visibility,dim(cartesian_profile)[2]*bin_width-incoming_distance)
        homogeneous_slant_visibility <- min(homogeneous_slant_visibility,dim(cartesian_profile)[2]*bin_width-incoming_distance)
      }
    }
    if (verbose | (!verbose && !output_files))
    {      
      if (incoming)
      {
        if (horizontal_visibility < incoming_distance)
        {
          message("Horizontal incoming visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m: ",horizontal_visibility," m.")
        } else {
          message("Incoming object at a height of ",incoming_height," m and distance of ",incoming_distance," m has optical contact with overhead airspace.")
        }
      } else {
        message("Horizontal visibility of outcoming object at a height of ",incoming_height," m and distance of ",incoming_distance," m: ",horizontal_visibility," m.")
      }
      if (vertical_visibility[1] < incoming_height)
      {
        message("Vertical visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m: ",vertical_visibility[1]," m.")
        message(median_slant_visibility)
        message(minimum_slant_visibility)
        message(progressive_slant_visibility)
        message(homogeneous_slant_visibility)
      } else {
        message(c("Outcoming","Incoming")[as.integer(incoming)+1]," object at a height of ",incoming_height," m and distance of ",incoming_distance," m has optical contact with ground.")
        message("Median slant visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m: ",median_slant_visibility," m.")
        message("Minimum slant visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m: ",minimum_slant_visibility," m.")
        message("Progressive slant visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m: ",progressive_slant_visibility," m.")
        message("Homogeneous slant visibility from a height of ",incoming_height," m and distance of ",incoming_distance," m: ",homogeneous_slant_visibility," m.")
      }
    }
  }
  
  #if output files are requested, two tables are generated: one containing the cartesian extinction profile and another containing the visibility ranges of the aircraft
  if (output_files)
  {
    if (is.null(model))
      model <- "no_model"
    message("Writing output files to disk...")
    if (!file.exists(paste("Cartesian_extinction_profile_",model,".txt.gz",sep="")))
    {
      message(paste("Cartesian_extinction_profile_",model,sep="")," not present. Please wait for first-time creation and compression.")
      write.table(cartesian_profile,file=paste("Cartesian_extinction_profile_",model,".txt",sep=""),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
      system(paste("gzip ","Cartesian_extinction_profile_",model,".txt",sep=""))      
    }
    if (!is.null(incoming_height) && !is.null(incoming_distance))
    {
      if (incoming)
      {
        write.table(c(horizontal_visibility,vertical_visibility[1],median_slant_visibility,minimum_slant_visibility,progressive_slant_visibility,homogeneous_slant_visibility),file=paste("Incoming_object_visibility_",model,".txt",sep=""),quote=FALSE,sep="\t",col.names=paste("Incoming height: ",incoming_height," m, incoming distance: ",incoming_distance," m.",sep=""),row.names=c("Horizontal visibility: ","Vertical visibility: ","Median slant visibility:","Minimum slant visibility:","Progressive slant visibility:","Homogeneous slant visibility:"))
      } else {
        write.table(c(horizontal_visibility,vertical_visibility[1],median_slant_visibility,minimum_slant_visibility,progressive_slant_visibility,homogeneous_slant_visibility),file=paste("Outcoming_object_visibility_",model,".txt",sep=""),quote=FALSE,sep="\t",col.names=paste("Outcoming height: ",incoming_height," m, outcoming distance: ",incoming_distance," m.",sep=""),row.names=c("Horizontal visibility: ","Vertical visibility: ","Median slant visibility:","Minimum slant visibility:","Progressive slant visibility:","Homogeneous slant visibility:"))
      }
    }
  }
  #return(cartesian_profile)
}
