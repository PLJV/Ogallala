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
#' using ofr98-393, make a contour line to raster product
#' @export
generateBaseElevationRaster <- function(s=NULL){
  # contourRaster <- raster(xmn=contourBBox[1,1],
  #                       xmx=ceiling(contourBBox[1,2]),
  #                       ymn=contourBBox[2,1],
  #                       ymx=ceiling(contourBBox[2,2]),
  #                       crs = CRS("+init=epsg:26918"),
  #                       resolution = c(12,12))
  # grid_pts <- as(s, 'SpatialPointsDataFrame')
  # g <- gstat::gstat(id="ELEV", formula = ELEV~1, data=p, nmax=7, set=list(idp = .5))
  # interpDEM <- gstat::interpolate(contourRaster, g) # idw interpolation
  return(NA)
}
#' calculate saturated thickness for a series of well points and the base elevation of
#' the aquifer (as returned by ogallala::generateBaseElevationRaster())
#' @export
calculateSaturatedThickness <- function(wellPts=NULL,baseRaster=NULL){

}
