#
# Workflow that generates a basic aquifer saturated thickness raster
# for a single year of well depth data
#
# Author : kyle.taylor@pljv.org
#

if(!require(Ogallala)) {
  install.packages("devtools", repos="https://cran.revolutionanalytics.com")
  stopifnot(require('devtools'))
  devtools::install_github("PLJV/Ogallala")
  stopifnot(require("Ogallala"))
}
if(!require('raster')){
  install.packages("raster", repos="https://cran.revolutionanalytics.com")
  stopifnot(require('raster'))
}
if(!require('rgdal')){
  install.packages("rgdal", repos="https://cran.revolutionanalytics.com")
  stopifnot(require('rgdal'))
}

# 2013 is the most recent year that the USGS has published. See:
# https://ne.water.usgs.gov/ogw/hpwlms/data.html
year <- 2013 

# define our regional boundary as the delineation of the aquifer
boundary <- get(
    'boundary', 
    envir=asNamespace('Ogallala')
  )

# fetch our input well data
wellPoints <- Ogallala:::unpackWellPointData(Ogallala:::scrapeWellPointData(
  years=year
  ))

# let's use the included base elevation and surface_elevation by default
# we could swap these out for our own if we wanted
base_elevation <- get(
    'base_elevation',
    envir=asNamespace('Ogallala')
  )

surface_elevation <- get(
    'surface_elevation',
    envir=asNamespace('Ogallala')
  )

# calculate saturated thickness field and adds it to our well points
# dataset
wellPoints <- Ogallala::calc_saturated_thickness(
    wellPts=wellPoints,
    baseRaster=base_elevation,
    surfaceRaster=surface_elevation,
    convert_to_imperial=T
  )

wellPoints <- wellPoints[!is.na(wellPoints@data$saturated_thickness),]

# polynomial trend surface from lat + lon + surface elev + bedrock elev
predictorData <- raster::stack(surface_elevation, base_elevation)
  names(predictorData) <- c("surf_elev","base_elev")

# create KNN smoothed field (literally averages sat_thickness values
# for each point using that point's neighbors -- from V. McGuire.

wellPoints <- Ogallala:::knn_point_smoother(
    wellPoints, 
    field="saturated_thickness",
    k=3 # explore K=2,3,4,etc...
  )

# train our interpolators

polynomialTrendSurface <- Ogallala::calc_polynomial_trend_surface(
    wellPoints,
    predRaster=predictorData, 
    field="saturated_thickness_smoothed",
    order=3    
  )

# check the fit of our model -- do our polynomial terms have a large effect?
# does the map look like a candy-cane? If so, reduce the order of our polynomials
# until the results look reasonable

summary(polynomialTrendSurface$m)

# are we predicting outside of the theoretical max sat. 
# thickness of the aquifer?
if(
  cellStats(polynomialTrendSurface$raster,'max') >
  cellStats(Ogallala:::meters_to_feet(surface_elevation-base_elevation), 'max')
){
  warning(
    "some of our predictions are outside of the theoretical",
    " max thickness thickness of the aquifer"
  )
}  
if(cellStats(polynomialTrendSurface$raster,'max') > 1300){
  warning(
    "max prediction is greater than the max reported depth of the aquifer",
    ". re-scaling using min-max normalization"
  )
  # re-shape our predictive interval so it is 0:1300 ft (from V McGuire)
  polynomialTrendSurface$raster <- Ogallala:::min_max_normalize(
    polynomialTrendSurface$raster)
}
# save to disk in the current working directory
cat(" -- writing raster to local directory\n")
writeRaster(
    polynomialTrendSurface$raster,
    paste(
      "saturated_thickness_",year,"_polynomial_trend.tif", 
      sep=""
    ),
    overwrite=T
  )

