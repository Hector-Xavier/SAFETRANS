arguments <- commandArgs(TRUE)
channels <- as.integer(strsplit(as.character(arguments[1]),split="_")[[1]])
is_scan <- as.logical(arguments[2])
scan_type <- as.character(arguments[3])
if (sum(scan_type==c("azimuth","elevation"))==0)
  stop("Incorrect scan type designation. Correct designations are: azimuth, elevation")
if (arguments[4]!="null")
{
  model <- as.character(arguments[4])
} else {
  model <- NULL
}
if (!is.null(model))
  if (sum(model==c("urban-rural","maritime","angstrom_exponent","all"))==0)
    stop("Incorrect model designation. Correct designations are: urban-rural, maritime, angstrom_exponent, null, all")
wavelength <- as.integer(arguments[5])
if (sum(wavelength==c(355,1064,1550))==0)
  stop("Incorrect lidar wavelength designation. Supported wavelengths (in nm): 355, 1064, 1550")
if (sum(model==c("urban-rural","maritime"))==1 && wavelength==355)
    stop("Empirical conversion available only for infrared wavelengths. Acceptable model designations for 355 nm: angstrom_exponent, null, all")
incoming <- as.logical(arguments[6])
if (arguments[7]!="null")
{
  height <- as.numeric(arguments[7])
} else {
  height <- NULL
}
if (arguments[8]!="null")
{
  distance <- as.numeric(arguments[8])
} else {
  distance <- NULL
}
output <- as.logical(arguments[9])
verbose <- as.logical(arguments[10])

source("scripta.R")
directory <- "/data"
setwd(directory)
if (output)
{
  if(!file.exists(paste(getwd(),"Output",sep="/")))
    dir.create(paste(getwd(),"Output",sep="/"))
  setwd(paste(directory,"Output",sep="/"))
}

if (!is_scan)
{
  radial_visibility_profile(extinction_profile=scanning_profile_extinction(directory,channels,is_scan,1,verbose,output),is_scan,model,wavelength,output,verbose)
} else {
  message("Processing measurements as a set of ",scan_type," measurements.")
  if (file.exists(paste("Radial_extinction_coefficients",wavelength,"nm.txt",sep="_")))
  {
    message("Radial extinction coefficent profile present. Loading data...")
    extinction <- read.table(file=paste("Radial_extinction_coefficients",wavelength,"nm.txt",sep="_"),header=TRUE,quote="")
  } else {
    if (!verbose)
      message("Importing data and calculating the radial extinction coefficent profile.")
    extinction <- scanning_profile_extinction(directory,channels,is_scan,1,verbose,output)
  }
  if (scan_type=="azimuth")
  {
    if (!verbose)
      message("Calculating the radial visibility profile.")
    if (model=="all")
    {
      radial_visibility_profile(extinction,is_scan,NULL,wavelength,output,verbose)
      if (wavelength!=355)
      {
        radial_visibility_profile(extinction,is_scan,"maritime",wavelength,output,verbose)
        radial_visibility_profile(extinction,is_scan,"urban-rural",wavelength,output,verbose)
      }
      radial_visibility_profile(extinction,is_scan,"angstrom_exponent",wavelength,output,verbose)
    } else {
      radial_visibility_profile(extinction,is_scan,model,wavelength,output,verbose)
    }
  } else {
    message("Extending radial extinction profile to cartesian coordinate system.")
    if (!is.null(height) && !is.null(distance))
      message("Calculating visibility ranges of ",c("outcoming","incoming")[as.integer(incoming)+1]," object.")
    if (model=="all")
    {
      cartesian_visibility_profile(extinction,NULL,wavelength,incoming,distance,height,output,verbose)
      if (wavelength!=355)
      {
      cartesian_visibility_profile(extinction,"maritime",wavelength,incoming,distance,height,output,verbose)
      cartesian_visibility_profile(extinction,"urban-rural",wavelength,incoming,distance,height,output,verbose)
      }
      cartesian_visibility_profile(extinction,"angstrom_exponent",wavelength,incoming,distance,height,output,verbose)
    } else {
      cartesian_visibility_profile(extinction,model,wavelength,incoming,distance,height,output,verbose)
    }
  }
}
