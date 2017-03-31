#' hidden scraping interface that will attempt to fetch well-point water data
scrapeWellPointData <- function( base_url="https://ne.water.usgs.gov/projects/HPA/", years=1995:(as.numeric(format(Sys.time(), "%Y"))) ){
  hrefs <- rvest::html(paste(base_url,"data.html",sep=""))
  hrefs <- as.character(rvest::html_nodes(hrefs,"a"))
  hrefs <- hrefs[grepl(hrefs,pattern="WL_ALL_")]
  hrefs <- hrefs[grepl(hrefs,
                       pattern=paste(years, collapse="|")
                       )]
  hrefs <- strsplit(hrefs,split="\"")
  hrefs <- unlist(lapply(hrefs,FUN=function(x) x[[2]]))

  zip_names <-  unlist(lapply(strsplit(hrefs,split="/"), FUN=function(x) x[length(x)]))

  toGet <- paste("https://ne.water.usgs.gov/projects/HPA/",hrefs,sep="")

  nExistingZips <- length(list.files(pattern="^WL_ALL_.*.zip$"))

  if(nExistingZips<length(years)){
    if(nExistingZips > 0){
      toGet <- list.files(pattern="^WL_ALL_.*.zip$")
        years <- years[which(!grepl(zips,pattern=paste(toGet,collapse="|")))]
          toGet <- zips[!grepl(zips, pattern=paste(toGet,collapse="|"))]
    }
    for(i in 1:length(toGet)){ utils::download.file(toGet[i],destfile=zip_names[i],quiet=T); cat('.') }
  }; cat("\n");
  return(list.files(pattern="^WL_ALL_.*.zip$"))
}
#' scrape bedrock elevation data
#  see: https://water.usgs.gov/GIS/metadata/usgswrd/XML/ofr98-393_aqbase.xml
#  note: [ELEV] field : base of aquifer elevation is in feet, not meters
scrapeBedrockElevation <- function(x=NULL){
  hrefs <- rvest::html("https://water.usgs.gov/GIS/metadata/usgswrd/XML/ofr98-393_aqbase.xml")
  hrefs <- rvest::xml_nodes(hrefs,"networkr")
  hrefs <- hrefs[grepl(hrefs,pattern="e00")]
  hrefs <- unlist(lapply(strsplit(as.character(hrefs),split=">|<"),function(x) x[3]))
  download.file(hrefs,destfile="ofr98-393.e00.gz")
}
#' hidden unpack and process well point data, returning to user as CSV
unpack_wl_zip <- function(x){
  default_projection <- "+proj=longlat +datum=NAD83 +no_defs" # strange projection (proj will figure it out)
  f <- utils::unzip(x,list=T)[,1];
    utils::unzip(x);
  f <-read.csv(f,sep="\t",comment.char="#",stringsAsFactors=FALSE)
    f[f<=-999] <- NA # looming NA values in the source CSV have inconsistently large values. Treat as NA's
       f[f>=9999] <- NA
  f <- sp::SpatialPointsDataFrame(coords=data.frame(x=f$long_dd_NAD83,y=f$lat_dd_NAD83),data=f,proj4string=sp::CRS(default_projection))
  return(sp::spTransform(f,sp::CRS(raster::projection("+init=epsg:2163"))))
}
#' unpack bedrock elevation data
unpack_bedrock <- function(x="ofr98-393.e00"){
  x = list.files(".",pattern=x)[1]
  if(!file.exists(x)) stop("couldn't find bedrock elevation contour data in cwd. use a better x= argument")
  if(grepl(x,pattern="gz")) R.utils::gunzip(x,overwrite=T)
  x <- rgdal::readOGR("ofr98-393.e00",verbose=F)
    return(sp::spTransform(x,sp::CRS(raster::projection("+init=epsg:2163"))))
}
