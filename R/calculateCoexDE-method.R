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
    noHKFlags     <- flagNotFullyExpressedGenes(objCOTAN)
    muEstimator   <- calculateMu(objCOTAN)
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
      set    <- set[!(set %in% getFullyExpressedGenes(objCOTAN))]
      complement <- setdiff(getGenes(objCOTAN), set) # genes not in this set
      complement <- complement[!(complement %in% getFullyExpressedGenes(objCOTAN))]

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
