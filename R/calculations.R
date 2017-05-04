#' hidden function for feet-> meters
feetToMeters <- function(x) x*0.3048
#' hidden function for meters -> feet
metersToFeet <- function(x) x/0.3048
#' generate a universal grid at the target resolution (500m) originally specified
#' by V. McGuire and others for their interpolation.
generateTargetRasterGrid <- function(s=NULL) {
  extent <- raster::extent(sp::spTransform(s,CRSobj=sp::CRS("+init=epsg:2163")))
  grid <- raster::raster(resolution=500,
                         ext=extent,
                         crs=sp::CRS("+init=epsg:2163"))
  return(grid)
}
#' using ofr98-393, make a contour line to raster product. The units on
#' bedrock elevation are feet (imperial). We are going to change them to
#' meters so they are consistent with surface DEMs from NED.
#' @param two_pass Boolean. Should we do a second pass of IDW interpolation to remove artifacts in
#' bedrock elevation.
#' @param feet_to_meters Boolean. Should we convert feet into meters?
#' @export
generateBaseElevationRaster <- function(s=NULL,targetRasterGrid=NULL,
                                        two_pass=T, feet_to_meters=T,
                                        mask=F){
  # sanity-check
  names(s) <- toupper(names(s))
  if(sum(grepl(names(s),pattern="ELEV"))==0) stop("no ELEV field found in the s= base countour shapefile provided.")
  if(is.null(targetRasterGrid)){
    cat(" -- using extent data from input Spatial* data to generate a 500m target grid\n")
    targetRasterGrid <- Ogallala:::generateTargetRasterGrid(s=s)
  }
  # local functions
  idw_interpolator <- function(pts,field="ELEV"){
    g <- gstat::gstat(id=field,
                      formula = as.formula(paste(field,"~1",sep="")),
                      data=pts)
    return(raster::interpolate(targetRasterGrid, g))
  }
  cat(" -- interpolating: ")
  cat(" -- pass one:")
  grid_pts <- as(s, 'SpatialPointsDataFrame')
  if(feet_to_meters){
    grid_pts$ELEV <- feetToMeters(grid_pts$ELEV)
  }
  bedrock <- idw_interpolator(grid_pts)
  # IDW is meant to be used on scattered data. Our contour lines were definately
  # not scattered and this results in localized artifacts in our interpolation
  # we are going to re-sample the raster again and re-do our interpolation
  # to try and smooth over these artifacts
  if(two_pass){
    cat(" -- pass two (resampling to remove artifacts):")
    resampled <- raster::sampleRandom(bedrock,size=999999,sp=T)
      names(resampled) <- "ELEV"
    bedrock <- idw_interpolator(resampled)
  }
  if(mask){
    cat(" -- masking\n")
    if(!file.exists(file.path("boundaries/ds543.zip"))){
      scrapeHighPlainsAquiferBoundary()
    }
    boundary <- sp::spTransform(unpackHighPlainsAquiferBoundary(),
                                sp::CRS(raster::projection(bedrock)))
    bedrock <- raster::mask(bedrock, boundary)
  }
  return(bedrock)
}
#' calculate saturated thickness for a series of well points and the base elevation of
#' the aquifer (as returned by ogallala::generateBaseElevationRaster())
#' @export
calculateSaturatedThickness <- function(wellPts=NULL,baseRaster=NULL){
  # sanity-check our input
  if(is.null(wellPts)) stop("wellPts= needs to be a SpatialPointsDataFrame specifying
                             well point data, as returned by unpackWellPointData()")
  if(is.null(baseRaster)){
    baseRaster = ogallala::generateBaseElevationRaster(s=wellPts)
  } else if(!inherits(baseRaster,'raster')){
    stop("baseRaster= argument must be a raster object, as returned by generateBaseElevationRaster()")
  }
  # calculate saturated thickness
  return(NA)
}
