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
  if (sum(model==c("urban-rural","maritime"))==0)
    stop("Incorrect model designation. Correct designations are: urban-rural, maritime, null")
incoming <- as.logical(arguments[5])
if (arguments[6]!="null")
{
  height <- as.numeric(arguments[6])
} else {
  height <- NULL
}
if (arguments[7]!="null")
{
  distance <- as.numeric(arguments[7])
} else {
  distance <- NULL
}
output <- as.logical(arguments[8])
verbose <- as.logical(arguments[9])

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
  radial_visibility_profile(extinction_profile=scanning_profile_extinction(directory,channels,is_scan,1,verbose,output),is_scan,model,output,verbose)
} else {
  message("Processing measurements as a set of ",scan_type," measurements.")
  if (file.exists("Radial_extinction_coefficients_1064_nm.txt"))
  {
    message("Radial extinction coefficent profile present. Loading data...")
    extinction <- read.table(file="Radial_extinction_coefficients_1064_nm.txt",header=TRUE,quote="")
  } else {
    if (!verbose)
      message("Importing data and calculating the radial extinction coefficent profile.")
    extinction <- scanning_profile_extinction(directory,channels,is_scan,1,verbose,output)
  }
  if (scan_type=="azimuth")
  {
    if (!verbose)
      message("Calculating the radial visibility profile.")
    radial_visibility_profile(extinction,is_scan,model,output,verbose)
  } else {
    message("Extending radial extinction profile to cartesian coordinate system.")
    if (!is.null(height) && !is.null(distance))
      message("Calculating visibility ranges of ",c("outcoming","incoming")[as.integer(incoming)+1]," object.")
    cartesian_visibility_profile(extinction,model,incoming,distance,height,output,verbose)
  }
}
