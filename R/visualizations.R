#' plot the z-standardized residual error for a model run
#' @export
plotResiduals <- function(pts, field='residls'){
  require(ggplot2)
  require(RColorBrewer)
  pts$lon <- pts@coords[,1]
  pts$lat <- pts@coords[,2]
  # compress the distribution so that outliers really stand out
  pts$res_plt <- quantileNormalize(pts@data[,field])
  ggplot(pts@data, aes(x=lon, y=lat, color=res_plt)) +
    geom_point(size=1.5, alpha=0.95) +
    geom_point(shape = 1,size = 1.5, color = "white", alpha=0.5) +
      xlab("longitude")+ ylab("latitude") + labs(color="Residuals (Z-score)") +
        scale_color_gradient(low="Red", high="Green")
}
