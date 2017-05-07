#' hidden function for feet-> meters [goo.gl/fguOPH]
feetToMeters <- function(x) x*0.3048
#' hidden function for meters -> feet [goo.gl/EmwvV1]
metersToFeet <- function(x) x*3.28084
#' generate a universal grid at the target resolution (500m) originally specified
#' by V. McGuire and others for their interpolation.
generateTargetRasterGrid <- function(s=NULL) {
  extent <- raster::extent(sp::spTransform(s,CRSobj=sp::CRS("+init=epsg:2163")))
  grid <- raster::raster(resolution=500,
                         ext=extent,
                         crs=sp::CRS("+init=epsg:2163"))
  return(grid)
}
#' hidden shortcut for gstat's idw interpolator
idw_interpolator <- function(pts,targetRasterGrid=NULL,field="ELEV"){
  pts <- pts[!is.na(pts@data[,field]),] # drop any NA values
  g <- gstat::gstat(id=field,
                    formula = as.formula(paste(field,"~1",sep="")),
                    data=pts)
  return(raster::interpolate(targetRasterGrid, g))
}
#' using ofr98-393, make a contour line to raster product. The units on
#' bedrock elevation are feet (imperial). We are going to change them to
#' meters so they are consistent with surface DEMs from NED.
#' @param two_pass Boolean. Should we do a second pass of IDW interpolation to remove artifacts in
#' bedrock elevation.
#' @param feet_to_meters Boolean. Should we convert feet into meters?
#' @export
generateBaseElevationRaster <- function(s=NULL,targetRasterGrid=NULL,
                                        two_pass=F, feet_to_meters=T,
                                        mask=F){
  # sanity-check
  names(s) <- toupper(names(s))
  if(sum(grepl(names(s),pattern="ELEV"))==0) stop("no ELEV field found in the s= base countour shapefile provided.")
  if(is.null(targetRasterGrid)){
    cat(" -- using extent data from input Spatial* data to generate a 500m target grid\n")
    targetRasterGrid <- Ogallala:::generateTargetRasterGrid(s=s)
      values(targetRasterGrid) <- 0
  }
  # build a K=5 NN index of random points across our target area
  # and populate with values informed by our base contour delineation
  cat(" -- interpolating\n")
  cat(" -- pass one: ")
  grid_pts <- as(s, 'SpatialPointsDataFrame')
    colnames(grid_pts@data) <- "ELEV"
      grid_pts <- raster::rasterize(grid_pts, targetRasterGrid, field="ELEV")
  base <- raster::resample(grid_pts,targetRasterGrid,method='bilinear',progress='text')
  if(feet_to_meters){
    base <- round(Ogallala:::feetToMeters(base),2)
  }
  if(mask){
    cat(" -- masking\n")
    if(!file.exists(file.path("boundaries/ds543.zip"))){
      scrapeHighPlainsAquiferBoundary()
    }
    boundary <- sp::spTransform(unpackHighPlainsAquiferBoundary(),
                                sp::CRS(raster::projection(bedrock)))
    base <- raster::mask(base, boundary)
  }
  return(base)
}
# use a KNN classifier fit to lat/lon to select a neighborhood of points around
# each well and calculates a summary statistic of your choosing (e.g., mean)
knnPointSmoother <- function(pts=NULL, field=NULL, k=4,fun=mean){
  index <- cbind(1:nrow(pts),
             FNN::get.knn(
               cbind(pts@data[,'base_elevation'], pts@coords), k=k)$nn.index)
  pts@data[,paste(field,"_smoothed",sep="")] <-
    apply(MARGIN=1,matrix(pts@data[as.vector(index),field],ncol=k+1),
          FUN=mean, na.rm=T)
  return(pts)
}
polynomialTrendSurface <- function(pts, order=4,
                                   field=NULL, predRaster=NULL){
  t <- cbind(pts@data[,field], pts@coords, pts$surface_elevation, pts$base_elevation)
    colnames(t) <- c(field,"longitude","latitude","surf_elev","base_elev")
      t <- data.frame(t)
  cat(" -- building a polynomial trend model\n")
  covs <- paste("poly(",colnames(t)[2:ncol(t)],",", order, ")",sep="")
    covs <- paste(covs,collapse="+")
      formula <- as.formula(paste(field,"~",covs,collapse=""))
  m <- glm(formula,data=na.omit(t))
  # if the user provided a rasterStack for making predictions, let's use it.
  if(!is.null(predRaster)){
    if(sum(!colnames(t)[2:ncol(t)] %in% names(predRaster)) > 0){
      # if we are missing predictors, are they latitude and longitude?
      if(sum(! colnames(t)[2:ncol(t)] %in% c(names(predRaster),"longitude","latitude")) == 2){
        cat(" -- calculating latitude and longitude")
        predRaster$latitude  <- init(predRaster,"y")
        predRaster$longitude <- init(predRaster,"x")
      }
    }
    cat(" -- projecting across regional extent:\n")
    polynomial_trend <- raster::predict(predRaster,m,progress='text',type="response")
    return(list(m=m,raster=polynomial_trend_out))
  }
  # return the model by default
  return(m)
}
splitToTrainingTestingDatasets <- function(pts=NULL,split=0.2){
  rows <- 1:nrow(pts)
  out_sample  <- sample(rows, size=split*nrow(pts))
  in_sample   <- rows[!rows %in% out_sample]
  return(list(training=pts[in_sample,],testing=pts[out_sample,]))
}
#' calculate saturated thickness for a series of well points and the base elevation of
#' the aquifer (as returned by ogallala::generateBaseElevationRaster())
#' @export
calculateSaturatedThickness <- function(wellPts=NULL,baseRaster=NULL,
                                        surfaceRaster=NULL, convert_to_imperial=T){
  # sanity-check our input
  if(is.null(wellPts)) stop("wellPts= needs to be a SpatialPointsDataFrame specifying
                             well point data, as returned by unpackWellPointData()")
  if(is.null(baseRaster)){
    baseRaster = Ogallala::generateBaseElevationRaster(s=wellPts)
  } else if(!inherits(baseRaster,'Raster')){
    stop("baseRaster= argument must be a raster object, as returned by generateBaseElevationRaster()")
  }
  if(is.null(surfaceRaster)){
    surfaceRaster = Ogallala::scrapeNed(s=wellPts)
  } else if(!inherits(baseRaster,'Raster')){
    stop("surfaceRaster= argument must be a raster object")
  }
  if(raster::projection(surfaceRaster) != raster::projection(baseRaster)){
    cat(" -- re-gridding surface raster to the resolution of our base raster\n")
    surfaceRaster <- projectRaster(surfaceRaster,
                                   to=baseRaster)
  }
  if(convert_to_imperial){
    # temporarily transform to metric so our depth observations
    # agree with our surface and base elevation datasets
    wellPts$well_depth_ft <- feetToMeters(wellPts$well_depth_ft)
    wellPts$lev_va_ft     <- feetToMeters(wellPts$lev_va_ft)
  }
  # extract surface and aquifer base information for our well points
  wellPts$surface_elevation <-
    raster::extract(x=surfaceRaster,
                    y=sp::spTransform(wellPts,CRS(projection(surfaceRaster))),
                    method='bilinear'
                   )
  wellPts$base_elevation <-
    raster::extract(x=baseRaster,
                    y=sp::spTransform(wellPts,CRS(projection(baseRaster))),
                    method='bilinear'
                   )
  # calculate saturated thickness
  depth_to_base <-
         (wellPts$surface_elevation - wellPts$base_elevation)
  wellPts$saturated_thickness <-
    depth_to_base - wellPts$lev_va_ft
  # units are always reported in (imperial) feet of saturated thickness
  wellPts$saturated_thickness <- if(using_metric) round(metersToFeet(wellPts$saturated_thickness),2) else round(wellPts$saturated_thickness,2)

  # some areas will have a non-sense "less than 0" depth_to_base
  # value -- assume that saturated thickness is zero in these places
  wellPts$saturated_thickness[wellPts$saturated_thickness<0] <- 0

  return(wellPts)
}
