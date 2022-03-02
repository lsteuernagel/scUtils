##########
### rasterize_ggplot
##########

#' Rasterize ggplot
#'
#' Makes use of scattermore::geom_scattermore to substitute GeomPoint layer in ggplot objects with an equivalent rasterized version.
#' Allows to manually modify the resolution of rasterized plots from Seurat plots!
#'
#' @param plot ggplot2 or patchwork object. For example output of Seurat's DimPlot or FeaturePlot
#' @param pixel_raster integer: number of pixels passed to pixels argument (x and optionall y) from scattermore::geom_scattermore. Defaults to 1024
#' @param pixel_raster_y pixel_raster_y: to use a different y pixels than x. Defaults to NULL which will use the value from pixel_raster in x and y
#' @param interpolate whether to linearly interpolate the image. Defaults to FALSE
#' @param pointsize Radius of rasterized point from scattermore. defaults to 1
#'
#' @return ggplot object with GeomPoint layers changed to GeomScattermore with rastrized version
#'
#' @export
#'
#' @import ggplot2 scattermore

rasterize_ggplot = function(plot,pixel_raster = 1024,pixel_raster_y = NULL,interpolate=FALSE,pointsize = 1){
  # if NULL will use pixel_raster else use different for y
  if(is.null(pixel_raster_y)){pixel_raster_y = pixel_raster}
  # get the geom point mapping
  layer_idx_all = which(sapply(plot$layers, function(x) class(x$geom)[1])=="GeomPoint")
  # check if any GeomPoint mapping exists
  if(length(layer_idx_all)==0){
    stop("Error: Cannot find points. 'rasterize_ggplot' expects to find at least one valid 'GeomPoint' layer to build GeomScattermore from.")
  }
  rasterized_plot = plot
  # in case there are multiple GeomPoint mappings: do for each:
  for(layer_idx in layer_idx_all){
    geom_point_layer = plot$layers[[layer_idx]]
    # make a plot with a rasterized geom
    rasterized_plot = rasterized_plot + scattermore::geom_scattermore(
      mapping = aes_string(
        x = as_label(geom_point_layer$mapping$x),
        y = as_label(geom_point_layer$mapping$y),
        color = paste0("`", as_label(geom_point_layer$mapping$colour), "`"),
        shape = as_label(geom_point_layer$mapping$shape),
        alpha = as_label(geom_point_layer$mapping$alpha)
      ),
      interpolate = interpolate,
      pointsize = pointsize,
      pixels = c(pixel_raster,pixel_raster_y)
    )
    # move into geom_point_mapping level in object (to preserve order of layers)
    rasterized_plot$layers[[layer_idx]] = rasterized_plot$layers[[length(rasterized_plot$layers)]]
    rasterized_plot$layers[[length(rasterized_plot$layers)]] = NULL
    #check if the new geom has non-default data that should be added
    if(length(geom_point_layer$data)>0){
      rasterized_plot$layers[[layer_idx]]$data = geom_point_layer$data
    }
  }
  return(rasterized_plot)
}
