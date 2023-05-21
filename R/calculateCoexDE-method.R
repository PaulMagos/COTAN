#' calculateCoexDE
#'
#' This function estimates the coex DE matrix '
#' @param objCOTAN A COTAN object
#' @param geneSets list of lists, every list is a set of genes
#'
#' @return a list with coex DE matrix and P-value
#'
#' @export
#'
#' @rdname calculateCoexDE
#'
setMethod(
  "calculateCoexDE",
  "COTAN",
  function(objCOTAN, geneSets){
    cat("data preparation..\n")

    # drop housekeeping genes
    noHKFlags     <- flagNotHousekeepingGenes(objCOTAN)
    muEstimator   <- estimateMu(objCOTAN)
    muEstimator   <- muEstimator[noHKFlags, ]
    zeroOneMatrix <- getZeroOneProj(objCOTAN)
    zeroOneMatrix <- zeroOneMatrix[noHKFlags, ]
    cells         <- getCells(objCOTAN)

    # Matrix probablity of 0
    probOfZero <- funProbZero(getDispersion(objCOTAN), muEstimator)

    # setData and setPVal will be returned as result
    setData <- as.data.frame(matrix(nrow = getNumCells(objCOTAN), ncol = 0))
    rownames(setData) <- getCells(objCOTAN)
    setPVal <- as.data.frame(matrix(nrow = getNumCells(objCOTAN), ncol = 0))
    rownames(setPVal) <- getCells(objCOTAN)

    for(setName in names(geneSets)){
      cat("working on set ", setName, "\n")

      set    <- unlist(geneSets[setName]) # needed to work with the list
      set    <- set[!(set %in% getHousekeepingGenes(objCOTAN))]
      complement <- setdiff(getGenes(objCOTAN), set) # genes not in this set
      complement <- complement[!(complement %in% getHousekeepingGenes(objCOTAN))]

      stopifnot("ERROR. Some gene in the set are not present in the COTAN object!"
                <- all(set %in% getGenes(objCOTAN)))

      # observed fields of contingency tables
      observedSetN    <- as.matrix(colSums(zeroOneMatrix[set, ] == 0))
      observedSetY    <- as.matrix(colSums(zeroOneMatrix[set, ] == 1) )
      observedComplementN <- as.matrix(colSums(zeroOneMatrix[complement, ] == 0))
      observedComplementY <- as.matrix(colSums(zeroOneMatrix[complement, ] == 1))

      # estimate fields of contingency tables
      extimateSetN    <- as.matrix(colSums(probOfZero[set, ]))
      extimateSetY    <- as.matrix(colSums(1 - probOfZero[set, ]))
      extimateComplementN <- as.matrix(colSums(probOfZero[complement, ]))
      extimateComplementY <- as.matrix(colSums(1 - probOfZero[complement, ]))

      # (in) set value 0 (No)
      inNo <- (observedSetN - extimateSetN) ** 2 /
        pmax(1, extimateSetN)

      # (out) complement value 0 (No)
      outNo <- (observedComplementN - extimateComplementN) ** 2 /
        pmax(1, extimateComplementN)

      # (in) set value 0 (Yes)
      inYes <- (observedSetY - extimateSetY) ** 2 /
        pmax(1, extimateSetY)

      # (out) complement value 0 (yes)
      outYes <- (observedComplementY - extimateComplementY) ** 2 /
        pmax(1, extimateComplementY)

      S <- inNo + outNo + inYes + outYes

      # clean memory. Now only S (the sum) is useful
      rm(inNo, outNo, inYes, outYes)
      gc()

      if (any(is.na(S))) {
        print(paste("Error: some NA in matrix S",
                    which(is.na(S), arr.ind = T)))
        break()
      }

      # P-value
      PValue <- as.data.frame(pchisq(as.matrix(S), df <- 1, lower.tail <- F))

      # COEX calculation
      coex <- (observedSetY - extimateSetY) / extimateSetY +
        (observedComplementN - extimateComplementN) / extimateComplementN -
        (observedComplementY - extimateComplementY) / extimateComplementY -
        (observedSetN - extimateSetN) / extimateSetN

      sumForDiv <- (1 / pmax(1, extimateComplementY) +
                      1 / pmax(1, extimateComplementN) +
                      1 / pmax(1, extimateSetY) +
                      1 / pmax(1, extimateSetN))

      # clean memory. Now only coex and sumForDivision are useful
      rm(observedComplementY, observedComplementN, observedSetY, observedSetN)
      rm(extimateComplementY, extimateComplementN, extimateSetY, extimateSetN)
      gc()

      coex <- coex / sqrt(sumForDiv)
      coex <- coex/sqrt(getNumCells(objCOTAN))
      coex <- as.data.frame(coex)

      colnames(coex)    <- paste0("set.", setName)
      colnames(PValue)  <- paste0("set.", setName)
      setData <- cbind(setData, coex)
      setPVal <- cbind(setPVal, PValue)

      gc()
    }

    return(list("sets" = setData, "pval" = setPVal))
  }
)

#' findSetsPoints
#'
#' this function subdivides the points of 'MatrixPoints' (arranged on the rows),
#' into sets based on the Euclidean distance that divides them.
#' @param matrixPonts:     matrix of points arranged in rows
#' @param maxDist:         maximum distance between two points in a set
#' @param minNumPointsSet: minimum number of points in a set
#' @param maxNumPointsSet: maximum number of points in a set
#' @param permutation:     if true, a permutation of the points is performed
#'
#' @export
#'
#' @return setsPoints:      output. list of lists, the names of the list are
#' the centers identified. The elements of one of the lists are the points
#' around the center of that list
#'
findSetsPoints <- function(matrixPoints,
                           maxDist = 15,
                           minNumPointsSet = 25,
                           maxNumPointsSet = 100,
                           permutation = FALSE){
  setsPoints <- vector(mode = "list")

  if(permutation){
    possibleCenters <- sample(rownames(matrixPoints))
  } else {
    possibleCenters <- rownames(matrixPoints)
  }

  repeat{
    center <- possibleCenters[1]

    inTheSet <- rep(Inf, nrow(matrixPoints))
    names(inTheSet) <- rownames(matrixPoints)

    cat("Working on center ", center, "\n")

    for(point in rownames(matrixPoints)){
      # Euclidean distance
      distance <- sqrt(sum((matrixPoints[center, ] - matrixPoints[point, ])^2))

      if(distance <= maxDist){
        inTheSet[point] <- distance
      }
    }

    numberNeighboringPoints <- sum(inTheSet != Inf)

    if(numberNeighboringPoints < minNumPointsSet){
      # the current center is not used as a center
      possibleCenters <- possibleCenters[2:length(possibleCenters)]
    } else {
      setSize         <- min(maxNumPointsSet, numberNeighboringPoints)
      # the closest points are inserted into the set
      namesOfClosest  <- names(sort(inTheSet))[1:setSize]
      possibleCenters <- setdiff(possibleCenters, namesOfClosest)

      setsPoints[[center]] <- namesOfClosest
    }

    if(is_empty(possibleCenters)){
      break
    }
  }

  return(setsPoints)
}


#' plotSetsOfPoints
#'
#' plot the overlap of the points in the sets. It is necessary that the
#' point names of all sets are in the names of 'x'.
#'
#' @param x vector of the first component of the points
#' @param y vector of the second component of the points
#' @param setsPoints list of lists. Each list has a name of 'x' as its name.
#' Each list contains a subset of the names of 'x'
#'
#' @returns ggplot
#'
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_gradient
#'

# @examples
# plotSetsOfPoints(pca[,1], pca[,2], SP)

plotSetsOfPoints <- function(x, y, setsPoints){
  namePoints         <- names(x)
  colorsValue        <- rep(0, length(namePoints))
  names(colorsValue) <- namePoints
  for(center in names(setsPoints)){
    colorsValue[setsPoints[[center]]] <- colorsValue[setsPoints[[center]]] + 1
  }
  colorsValue[names(setsPoints)] <- max(colorsValue) + 1

  plotDF <- data.frame(x = x, y = y, color = colorsValue)
  myPlot <- ggplot(plotDF, aes(x = x,
                               y = y,
                               color = color)) +
    geom_point() +
    scale_colour_gradient(low="grey100", high="red")

  return(list("plot" = myPlot, "overlapping" = colorsValue))
}


