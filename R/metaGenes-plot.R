
#' @details `metaGenesPlot()` shows the overlap of the points in the
#'   *meta-genes* in the given \eqn{(x, y)} space (usually the PCA)
#'
#'
#' @param x an `array` giving the first coordinate of the points
#' @param y an `array` giving the second coordinate of the points
#' @param metaGenes a `list` of *meta-genes*, such as the result of
#'   [defineMetaGenes()]. The union of all names in `metaGenes` must be a subset
#'   of the names in the coordinates
#'
#' @returns `metaGenesPlot()` returns a `list` with:
#'  * "plot" the wanted `ggplot`
#'  * "overlapping" an `array` that for each point gives how many
#'    *meta-genes* it belongs to; centers are returned with the highest value
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_gradient
#'
#' @importFrom rlang set_names
#'
#' @export
#'
#' @rdname MetaGenes

# @examples
# metaGenesPlot(pca[,1], pca[,2], SP)

metaGenesPlot <- function(x, y, metaGenes) {
  pointNames <- names(x)
  colorsValue <- set_names(rep(0, length(pointNames)), pointNames)

  for (center in names(metaGenes)) {
    metaGene <- metaGenes[[center]]
    colorsValue[metaGene] <- colorsValue[metaGene] + 1
  }
  colorsValue[names(metaGenes)] <- max(colorsValue) + 1

  plotDF <- data.frame(x = x, y = y, color = colorsValue)
  myPlot <- ggplot(plotDF) +
    geom_point(aes(x = x, y = y, color = color, size = 0.5)) +
    scale_colour_gradient(low = "grey100", high = "red")

  return(list("plot" = myPlot, "overlapping" = colorsValue))
}
