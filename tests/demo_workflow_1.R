#
# Workflow 2009 Aquifer saturated thickness validation and ensemble prediction
# generator
#
# Author : kyle.taylor@pljv.org
#

require(Ogallala)

# define our regional boundary
boundary <- Ogallala:::scrapeHighPlainsAquiferBoundary()
  boundary <- Ogallala:::unpackHighPlainsAquiferBoundary(boundary)

# fetch our input well data
wellPoints <- Ogallala:::scrapeWellPointData(years=2009)
  wellPoints <- Ogallala:::unpackWellPointData(wellPoints)

# build an aquifer base raster
base_elevation <- Ogallala:::scrapeBaseElevation()
  base_elevation <- Ogallala:::unpackBaseElevation()
    base_elevation <- Ogallala:::generateBaseElevationRaster(base_elevation)
writeRaster(base_elevation, "base_elevation.tif")

# calculate saturated thickness
wellPts <- Ogallala:::calculateSaturatedThickness(wellPts=wellPts,
                                       baseRaster=base_elevation,
                                       convert_to_imperial=T)

# create KNN smoothed field
wellPts <- Ogallala:::knnPointSmoother(wellPts, field="saturated_thickness")

# train our interpolators

# topogrid was pre-calculated in ArcGIS (below)

# idw
inverse_distance <- Ogallala:::idw_interpolator(wellPts,
                      targetRasterGrid=base_elevation,
                        field="saturated_thickness")

writeRaster(inverse_distance, "saturated_thickness_09_idw.tif")

inverse_distance_w_knn <- Ogallala:::idw_interpolator(wellPts,
                            targetRasterGrid=base_elevation,
                              field="saturated_thickness_smoothed")

writeRaster(inverse_distance_w_knn, "saturated_thickness_09_knn_idw.tif")

# polynomial trend
polynomial_trend <- Ogallala:::polynomialTrendSurface(wellPts, field="saturated_thickness_smoothed")
  writeRaster(polynomial_trend$raster, "saturated_thickness_09_polynomial_trend.tif")

# do some cross-validation

inverse_distance       <- raster("saturated_thickness_09_idw.tif")
inverse_distance_w_knn <- raster("saturated_thickness_09_knn_idw.tif")
polynomial_trend       <- raster("saturated_thickness_09_polynomial_trend.tif")
topogrid               <- projectRaster(raster("satThick.2009.vMcguire.tif"),
                            to=inverse_distance_w_knn)

# build two predictive ensembles to test, one with topogrid and one with polynomial trend
ensemble_tp <- stackApply(raster::stack(inverse_distance_w_knn,topogrid), fun=median, indices=1)
ensemble_pt <- stackApply(raster::stack(inverse_distance_w_knn,polynomial_trend), fun=median, indices=1)

# K-fold validation using random well points subsampled from the
# dataset. IDW (without neighbor weighting) will be our null model. This is
# a know-nothing interpolation with the original saturated thickness estimate
# calculated for each well. The other models being tested (including the ensembles)
# contain increasingly more information.
# Lastly, this validation testing will be slightly biased, because we aren't actually
# re-training models (topogrid isn't easily scriptable for this task). This might slightly
# inflate the performance of IDW algorithms, which we might expect to perform more poorly
# with fewer points. Regardless, all models will be biased the same way (i.e., have seen the
# full 2009 well dataset).

k = 100
overall <- data.frame(matrix(NA,nrow=k,ncol=6))
  colnames(overall) <- c("topogrid","idw","idw_w_knn",
                           "polynomial","ensemble_pt",
                             "ensemble_tp")

residuals <- data.frame(matrix(NA,nrow=nrow(pts)*0.2,ncol=4))

topogrid_raw_residual_error         <- vector()
polynomial_trend_raw_residual_error <- vector()
idw_w_knn_raw_residual_error        <- vector()

cat(" -- resampling and validating (pseudo-k-folds validation): ")
for(i in 1:k){
  run <- splitToTrainingTestingDatasets(pts)
  tg <-
    raster::extract(topogrid,
      spTransform(run$testing,CRS(projection(topogrid))))
  idw <-
    raster::extract(inverse_distance,
      spTransform(run$testing,CRS(projection(inverse_distance))))
  idw_w_knn <-
    raster::extract(inverse_distance_w_knn,
      spTransform(run$testing,CRS(projection(inverse_distance_w_knn))))
  pt <-
    raster::extract(polynomial_trend,
      spTransform(run$testing,CRS(projection(polynomial_trend))))
  e_pt <-
    raster::extract(ensemble_pt,
      spTransform(run$testing,CRS(projection(ensemble_pt))))
  e_tp <-
    raster::extract(ensemble_tp,
      spTransform(run$testing,CRS(projection(ensemble_tp))))

  # we are going to keep a running log of our residual error with
  # topogrid to see if we can find a directional bias
  topogrid_raw_residual_error <-
    append(topogrid_raw_residual_error,
      run$testing$saturated_thickness-tg)

  polynomial_trend_raw_residual_error <-
    append(polynomial_trend_raw_residual_error,
      run$testing$saturated_thickness-pt)

  idw_w_knn_raw_residual_error <-
    append(idw_w_knn_raw_residual_error,
      run$testing$saturated_thickness-idw_w_knn)

  residuals[,1] <- (run$testing$saturated_thickness-tg)^2
    overall[i,1] <- mean(residuals[,1],na.rm=T)
  residuals[,2] <- (run$testing$saturated_thickness-idw)^2
    overall[i,2] <- mean(residuals[,2],na.rm=T)
  residuals[,3] <- (run$testing$saturated_thickness-idw_w_knn)^2
    overall[i,3] <- mean(residuals[,3],na.rm=T)
  residuals[,4] <- (run$testing$saturated_thickness-pt)^2
    overall[i,4] <- mean(residuals[,4],na.rm=T)
  residuals[,5] <- (run$testing$saturated_thickness-e_pt)^2
    overall[i,5] <- mean(residuals[,5],na.rm=T)
  residuals[,6] <- (run$testing$saturated_thickness-e_tp)^2
    overall[i,6] <- mean(residuals[,6],na.rm=T)

  cat(paste("[",i,"/",k,"]",sep=""))
};
cat("\n")
# let's merge residuals from our last testing dataset into
# a shapefile for mapping
colnames(residuals) <- c("topogrid","idw","idw_w_knn",
                         "polynomial","ensemble_pt",
                           "ensemble_tp")
  run$testing@data <- cbind(run$testing@data,residuals)
    rgdal::writeOGR(run$testing,".", "well_pts_testing_validation",
                    overwrite=T, driver="ESRI Shapefile")
#
# report the winning prediction
#
winner <- names(which(round(colMeans(overall)) == min(round(colMeans(overall)))))
cat(" -- the winning algorithm is:", winner,"\n")
print(round(colMeans(overall)))

#
# Now Make Our Plots
#

#
# Plot 1 : Regional Predictive Ambiguity in Models vs. IDW Ensembles
#

dev.new()
par(mfrow=c(2,2), mai=c(0.1,0.55,0.25,1.25))
plot(topogrid,main="topogrid (2009)",cex=0.6, yaxt='n', xaxt='n')
plot(spTransform(boundary,CRS(projection(topogrid))),add=T,border="grey")
plot(polynomial_trend,main="polynomial trend (2009)",cex=0.6, yaxt='n', xaxt='n')
plot(spTransform(boundary,CRS(projection(polynomial_trend))),add=T,border="grey")
plot(ensemble_tp,main="ensemble topogrid (2009)",cex=0.6, yaxt='n', xaxt='n')
plot(spTransform(boundary,CRS(projection(ensemble_pt))),add=T,border="grey")
plot(ensemble_pt,main="ensemble polynomial trend (2009)",cex=0.6, yaxt='n', xaxt='n')
plot(spTransform(boundary,CRS(projection(ensemble_pt))),add=T,border="grey")

#
# Plot 2 : Residual error plots
#

dev.new()
hist(na.omit(topogrid_raw_residual_error),main="",xlab="residual error",cex=0.8,breaks=75,ylab="")
  grid()
    abline(col="red",v=mean(topogrid_raw_residual_error,na.rm=T))
dev.new()
hist(na.omit(polynomial_trend_raw_residual_error),main="",xlab="residual error",cex=0.8,breaks=75,ylab="")
  grid()
    abline(col="red",v=mean(polynomial_trend_raw_residual_error,na.rm=T))
