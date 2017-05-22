require(parallel)

years <- 2000:2014

buildPolynomialTrendEnsemble <- function(y=NULL, write=FALSE){
  require(sp)
  require(raster)
  require(rgdal)
  require(gstat)

  target <- paste("saturated_thickness_",y,"_ensemble_idw_polynomial_trend.tif",sep="")

  if(file.exists(target)){
    return(raster::raster(target))
  }

  surface_elevation <- raster::raster("surface_elevation.tif")
     base_elevation <- raster::raster("base_elevation.tif")

   # define our regional boundary
   boundary <- Ogallala:::scrapeHighPlainsAquiferBoundary()
     boundary <- Ogallala:::unpackHighPlainsAquiferBoundary(boundary)

  wellPoints <- sp::spTransform(Ogallala:::unpackWellPointData(Ogallala:::scrapeWellPointData(years=y))
                  , sp::CRS(raster::projection(base_elevation)))

  wellPoints <- Ogallala:::calculateSaturatedThickness(wellPts=wellPoints,
    baseRaster=base_elevation,
    surfaceRaster=surface_elevation,
    convert_to_imperial=T)

  wellPoints <- wellPoints[!is.na(wellPoints@data$saturated_thickness),]

  # create KNN smoothed field
  wellPoints <- Ogallala:::knnPointSmoother(wellPoints, field="saturated_thickness")
  inverse_distance_w_knn <- Ogallala:::idw_interpolator(wellPoints,targetRasterGrid=base_elevation,field="saturated_thickness_smoothed")

  # polynomial trend
  predictor_data <- raster::stack(surface_elevation, base_elevation)
    names(predictor_data) <- c("surf_elev","base_elev")
  polynomial_trend <- Ogallala:::polynomialTrendSurface(wellPoints,
    predRaster=predictor_data, field="saturated_thickness_smoothed")

  # ensemble
  ensemble_pt <- stackApply(raster::stack(inverse_distance_w_knn,polynomial_trend$raster), fun=median, indices=1)
  ensemble_pt <- raster::mask(ensemble_pt,
    sp::spTransform(boundary,sp::CRS(raster::projection(ensemble_pt))))
  if(write){
    writeRaster(ensemble_pt,
                paste("saturated_thickness_",y,"_ensemble_idw_polynomial_trend.tif",sep=""),
                progress='text')
  }
  return(ensemble_pt)
}

cl <- makeCluster(6)
out <- parallel::parLapply(cl, as.list(years),fun=buildPolynomialTrendEnsemble,write=T)
