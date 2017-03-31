#' scraping interface that will attempt to fetch well-point water data
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
#'unpack and process well point data, returning to user as CSV
unpack_wl_zip <- function(x){
  f <- utils::unzip(x,list=T)[,1];
    utils::unzip(x);
  f <-read.csv(f,sep="\t",comment.char="#",stringsAsFactors=FALSE)
    f[f<=-999] <- NA # looming NA values in the source CSV have inconsistently large values. Treat as NA's
       f[f>=9999] <- NA
  return(f)
}
