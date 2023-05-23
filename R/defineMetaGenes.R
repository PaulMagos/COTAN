#' Meta-genes: defintion and usage
#'
#' @description A *meta-gene* is a small group of genes deemed similar according
#'   to some given distance. This grouping is intended to make less noisy and
#'   more robust any procedure that would try to cluster cells based on coex
#'   information
#'
#' @details `defineMetaGenes()` groups together *similar* genes into
#'   *meta-genes* according to the given distance. Usually the distance is based
#'   on some function of the PCA of the genes
#'
#' @param points a `matrix` where each point is presented as a row
#' @param maxDist a maximum distance between two points inside a single
#'   *meta-gene*
#' @param minSize: the minimum size of a *meta-gene*
#' @param maxSize: the maximum size of a *meta-gene*
#' @param distance the distance method to use. Default is `"euclidian"`
#' @param permute: a Boolean. When `FALSE` the points are used in the given
#'   order as potential meta-genes centers, otherwise a random a permutation is
#'   performed
#'
#' @return `defineMetaGenes()`returns a `list` of *meta-genes*, each named after
#'   its center gene. A *metag-gene* is an array of gene names, thus making this
#'   similar to a [ClusterList]
#'
#' @export
#'
#' @importFrom parallelDist parDist
#'
#' @rdname MetaGenes
#'
defineMetaGenes <- function(points,
                            maxDist = 15.0,
                            minSize = 25L,
                            maxSize = 100L,
                            distance = "cosine",
                            permutation = FALSE) {
  logThis("Defining meta-genes - START", logLevel = 2L)

  metaGenes <- vector(mode = "list")

  pointsDist <- as.matrix(parDist(points, method = distance,
                                  diag = TRUE, upper = TRUE))

  possibleCenters <- rownames(matrixPoints)
  if (isTRUE(permutation)) {
    possibleCenters <- sample(possibleCenters)
  }

  repeat{
    center <- possibleCenters[[1L]]

    logThis(paste0("Working on center: ", center), logLevel = 3L)

    distances <- sort(pointsDist[center, ], decreasing = FALSE)

    numNeighbors <- sum(distances <= maxDist)

    if (numNeighbors < minNum) {
      # the current center is not used as a center
      possibleCenters <- possibleCenters[2L:length(possibleCenters)]
    } else {
      # the closest points are inserted into the set
      namesOfClosest  <- names(head(distances, n = min(maxNum, numNeighbors)))
      possibleCenters <- setdiff(possibleCenters, namesOfClosest)

      metaGenes[[center]] <- namesOfClosest
    }

    if (is_empty(possibleCenters)) {
      # no points left to process
      break
    }
  }

  return(metaGenes)
}
