#' hidden function for feet-> meters
feetToMeters <- function(x) x*0.3048
#' generate a universal grid at the target resolution (500m) originally specified
#' by V. McGuire and others for their interpolation.
generateTargetRasterGrid <- function(s=NULL) {
  extent <- raster::extent(sp::spTransform(s,CRSobj=sp::CRS("+init=epsg:2163")))
  grid <- raster::raster(resolution=500,
                         ext=extent,
                         crs=sp::CRS("+init=epsg:2163"))
  return(grid)
}
#' using ofr98-393, make a contour line to raster product. The units on elevation are feet (imperial).
#' @export
generateBaseElevationRaster <- function(s=NULL,targetRasterGrid=NULL){
  # sanity-check
  names(s) <- toupper(names(s))
  if(sum(grepl(names(s),pattern="ELEV"))==0) stop("no ELEV field found in the s= base countour shapefile provided.")
  if(is.null(targetRasterGrid)){
    cat(" -- no target raster grid specified. Will generate one from s= elevation data.\n")
    targetRasterGrid <- generateTargetRasterGrid(s=s)
  }

  grid_pts <- as(s, 'SpatialPointsDataFrame')
         g <- gstat::gstat(id="ELEV", formula = ELEV~1, data=grid_pts)

  bedrock <- interpolate(targetRasterGrid, g)

  return(bedrock)
}
#' calculate saturated thickness for a series of well points and the base elevation of
#' the aquifer (as returned by ogallala::generateBaseElevationRaster())
#' @export
calculateSaturatedThickness <- function(wellPts=NULL,baseRaster=NULL){
  # sanity-check our input
  if(is.null(wellPts)) stop("wellPts= needs to be a SpatialPointsDataFrame specifying
                             well point data, as returned by unpack_wl_zip()")
  if(is.null(baseRaster)){
    baseRaster = ogallala::generateBaseElevationRaster(s=wellPts)
  } else if(!inherits(baseRaster,'raster')){
    stop("baseRaster= argument must be a raster object, as returned by generateBaseElevationRaster()")
  }
  # calculate saturated thickness
  wellPts$
}
