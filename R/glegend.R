#' Save ggplot2 Legend
#' 
#' Stores or extracts a ggplot2 legend
#'
#' @param a.gplot a plot made with ggplot2
#' @export
glegend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}