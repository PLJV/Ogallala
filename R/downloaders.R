#' subset by grep match
grep_by <- function(x, pattern=NULL){
  x[grepl(x, pattern=pattern)]
}
#' grab all zips by default and parse through what the user requested
parse_hrefs <- function(url=NULL){
  a <- xml2::read_html(url)
    return(as.character(rvest::html_nodes(a,"a")))
}
#' grab the full string within a set of parantheses
parse_url <- function(x=NULL,pattern=NULL){
  if(grepl(x, pattern="\"")){
    split <- unlist(lapply(strsplit(x,split="/"), FUN=function(x) x[length(x)]))

  } else if(grepl(x, pattern="<|>")){
    split <- unlist(lapply(strsplit(as.character(x),split=">|<"),function(x) x[3]))
  }
  unlist(lapply(split,
                FUN = function(x) x[which(grepl(x, pattern = pattern))]
                ))
}
#' scrape hp aquifer boundary
scrapeHighPlainsAquiferBoundary <- function(base_url="https://water.usgs.gov/GIS/metadata/usgswrd/XML/ds543.xml"){
  if(!dir.exists("boundaries")){
    dir.create("boundaries")
  }
  hrefs <- xml2::read_html(base_url)
  hrefs <- rvest::xml_nodes(hrefs,"networkr")
  hrefs <- grep_by(hrefs,pattern="ds543")
  hrefs <- parse_url(hrefs,pattern="zip$")

  zip_name <-  unlist(lapply(strsplit(hrefs,split="/"), FUN=function(x) x[length(x)]))
    zip_name <- file.path(paste("boundaries/",zip_name,sep=""))
  utils::download.file(hrefs[1],destfile=zip_name,quiet=T);
  return(zip_name)
}
#' scrape well-point water data
scrapeWellPointData <- function( base_url="https://ne.water.usgs.gov/projects/HPA/", years=1995:(as.numeric(format(Sys.time(), "%Y"))) ){
  hrefs <- parse_hrefs(paste(base_url,"data.html",sep=""))
  hrefs <- grep_by(hrefs, pattern="WL_ALL_")
  hrefs <- grep_by(hrefs, pattern=paste(years, collapse="|"))
  hrefs <- parse_url(hrefs,pattern="zip$")

  zip_names <-  unlist(lapply(strsplit(hrefs,split="/"), FUN=function(x) x[length(x)]))

  toGet <- paste("https://ne.water.usgs.gov/projects/HPA/",hrefs,sep="")
  existingZips <- list.files(pattern="^WL_ALL_.*.zip$")

  if(length(existingZips) < length(toGet)){
    for(i in 1:length(existingZips)){
      toGet <- toGet[which(!grepl(toGet,pattern=existingZips[i]))]
    }
    for(i in 1:length(toGet)){
      utils::download.file(toGet[i],destfile=zip_names[i],quiet=T); cat('.')
    }
  }
  cat("\n");
  return(list.files(pattern="^WL_ALL_.*.zip$"))
}
#' scrape bedrock elevation data
#  see: https://water.usgs.gov/GIS/metadata/usgswrd/XML/ofr98-393_aqbase.xml
#  note: [ELEV] field : base of aquifer elevation is in feet, not meters
scrapeBedrockElevation <- function(base_url="https://water.usgs.gov/GIS/metadata/usgswrd/XML/ofr98-393_aqbase.xml"){
  hrefs <- xml2::read_html(base_url)
  hrefs <- rvest::xml_nodes(hrefs,"networkr")
  hrefs <- hrefs[grepl(hrefs,pattern="e00")]
  hrefs <- unlist(lapply(strsplit(as.character(hrefs),split=">|<"),function(x) x[3]))
  download.file(hrefs,destfile="ofr98-393.e00.gz")
}
#' hidden unpack and process well point data, returning to user as CSV
unpackWellPointData <- function(x=NULL){
  default_projection <- "+proj=longlat +datum=NAD83 +no_defs" # strange projection (proj will figure it out)
  f <- utils::unzip(x,list=T)[,1]; # get a file listing
       utils::unzip(x)             # actually decompress the file
  f <-read.csv(f,sep="\t",comment.char="#",stringsAsFactors=FALSE)
    f[f<=-999] <- NA # looming NA values in the source CSV have inconsistently large values. Treat as NA's
    f[f>=9999] <- NA
      f <- f[!(is.na(f$lat_dd_NAD83) | is.na(f$lat_dd_NAD83)),]
  f <- sp::SpatialPointsDataFrame(coords=data.frame(x=f$long_dd_NAD83,y=f$lat_dd_NAD83),data=f,proj4string=sp::CRS(default_projection))
  return(sp::spTransform(f,sp::CRS(raster::projection("+init=epsg:2163"))))
}
#' unpack bedrock elevation data
unpackBedrockElevation <- function(x="ofr98-393.e00"){
  x = list.files(".",pattern=x)[1]
  if(!file.exists(x)) stop("couldn't find bedrock elevation contour data in cwd. use a better x= argument")
  if(grepl(x,pattern="gz")) R.utils::gunzip(x,overwrite=T)
  x <- rgdal::readOGR("ofr98-393.e00",verbose=F)
    return(sp::spTransform(x,sp::CRS(raster::projection("+init=epsg:2163"))))
}
#' unpack aquifer boundary data
unpackHighPlainsAquiferBoundary <- function(x="ds543.zip"){
  x <- file.path(list.files(".",pattern="ds543.zip",recursive=T))
  if(!file.exists(x)) stop("couldn't find high plains boundary contour data in cwd. use a better x= argument")
  f <- utils::unzip(x,list=T,)[1,1];        # get a file listing
       utils::unzip(x,exdir="boundaries")             # actually decompress the file

  f <- gsub(f,pattern="[.]...", replacement="") # get layer name
  return(rgdal::readOGR("boundaries",f,verbose=F))
}
