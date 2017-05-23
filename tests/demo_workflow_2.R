#
# Workflow 2008-2013 ensemble prediction. I did some parallelization here to
# speed up raster generation across the time series, do some residual error
# calculations for each year, and do a grid-cell level time-series regression.
# This is just demonstrative. Looking longer term, that we will instead
# fit our regressions to the individual well point data
#
# Author : kyle.taylor@pljv.org
#
require(parallel)

years <- 2008:2013

buildPolynomialTrendEnsembleRasters <- function(y=NULL, write=FALSE, calc_residuals=FALSE){
  require(sp)
  require(raster)
  require(rgdal)
  require(gstat)

  target <- paste("saturated_thickness_",y,"_ensemble_idw_polynomial_trend.tif",sep="")

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

  # fetch spatially-weighted pseudo-zero values from the ofr99-266 report
  wellPoints <- Ogallala:::generatePseudoZeros(wellPoints)

  # create KNN smoothed field
  wellPoints <- Ogallala:::knnPointSmoother(wellPoints, field="saturated_thickness")
  if(!file.exists(target)){
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
    # write to disk, if asked
    if(write){
      writeRaster(ensemble_pt,
                  paste("saturated_thickness_",y,"_ensemble_idw_polynomial_trend.tif",sep=""),
                  progress='text')
    }
  } else {
    ensemble_pt <- raster::raster(target)
  }
  if(calc_residuals){
    wellPoints$residuals <- wellPoints$saturated_thickness_smoothed - raster::extract(ensemble_pt,wellPoints,df=F,sp=F)
    writeOGR(wellPoints,".",gsub(target,pattern="[.]tif$",replacement=""),driver="ESRI Shapefile", overwrite=T)
  }

  return(ensemble_pt)
}

plotResiduals <- function(pts){
  require(ggplot2)
  require(RColorBrewer)
  pts$lon <- pts@coords[,1]
  pts$lat <- pts@coords[,2]
  # compress the distribution so that outliers really stand out
  div <- round(diff(quantile((pts$residls),na.rm=T,p=c(0.5,0.75))))
  pts$res_plt <- (pts$residls)/div
  ggplot(pts@data, aes(x=lon, y=lat, color=res_plt)) +
    geom_point(size=1.5, alpha=0.95) +
    geom_point(shape = 1,size = 1.5, color = "white", alpha=0.5) +
      xlab("longitude")+ ylab("latitude") +
        scale_color_gradient(low="Green", high="Red")
}

cl <- makeCluster(6)
sat_thickness_2008_2013 <- parallel::parLapply(cl, as.list(years),fun=buildPolynomialTrendEnsembleRasters,write=T,calc_residuals=T)

mean_sat_thickness_2008_2013 <- stackApply(raster::stack(sat_thickness_2008_2013),fun=mean,indices=1)
sd_sat_thickness_2008_2013 <- stackApply(raster::stack(sat_thickness_2008_2013),fun=sd,indices=1)

writeRaster(mean_sat_thickness_2008_2013,"mean_sat_thickness_2008_2013.tif",overwrite=T)
writeRaster(sd_sat_thickness,"sd_sat_thickness_2008_2013.tif",overwrite=T)

contours <- raster::rasterToContour(sat_thickness_2008_2013,levels=c(30,50,100,150,300,500))
  writeOGR(contours,".","contour_mean_sat_thickness_2008_2013",driver="ESRI Shapefile",overwrite=T)

lm_intercept <- calc(raster::stack(sat_thickness_2008_2013), fun = function(x) {
  if (all(is.na(x)))
    return(NA)
  else
    return(coef(lm(x ~ years))[1])
})

writeRaster(lm_intercept,"intercept_change_sat_thickness_2008_2013.tif", overwrite=T)

lm_slope <- calc(raster::stack(sat_thickness_2008_2013), fun = function(x) {
  if (all(is.na(x)))
    return(NA)
  else
    return(coef(lm(x ~ years))[2])
})

writeRaster(lm_slope,"slope_change_sat_thickness_2008_2013.tif",overwrite=T)
