#  talysExforMapping - Convert Between EXFOR and TALYS
#  Copyright (C) 2019  Georg Schnabel
#  
#  talysExforMapping is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  talysExforMapping is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>

#' Create Subentry Handler
#'
#' @return List with functions:
#'         \describe{
#'           \item{getName()}{returns signature of the general handler}
#'           \item{getInfo()}{returns additional information about handler}
#'           \item{canMap(subent)}{returns a bool indicating whether \code{subent} can be mapped}
#'           \item{needs(subent)}{returns data table with extended output specification, see ...}
#'           \item{map(subent)}{TODO}
#'         }
#' @export
#' @import data.table Matrix
#'
createSubentHandler <- function(subentHandlerList) {


  getName <- function() {
    "handler_general"
  }


  getInfo <- function() {
    ""
  }


  canMap <- function(subent, quiet=TRUE) {

    hnd <- getResponsibleHandler(subent)

    if (is.null(hnd) && !quiet) {
      scores <- sapply(subentHandlerList, function(curHandler)
        curHandler$canMap(subent, quiet=TRUE, ret.score=TRUE,
                          masterHandler = masterHandler))
      ranking <- order(scores, decreasing = TRUE)
      for (idx in ranking[seq_len(min(length(ranking),3))]) {
        cat("--- Handler: ", subentHandlerList[[idx]]$getName(), " ---\n")
        subentHandlerList[[idx]]$canMap(subent, quiet=FALSE,
                                        masterHandler = masterHandler)
      }
    }
    !is.null(hnd)
  }


  needs <- function(subent, rowidcs=NULL) {

    hnd <- getResponsibleHandler(subent)
    if (is.null(hnd))
      stop("no appropriate handler available")
    hnd$needs(subent, rowidcs,
              masterHandler = masterHandler)
  }


  map <- function(subent, extout, rowidcs=NULL) {

    hnd <- getResponsibleHandler(subent)
    if (is.null(hnd))
      stop("no appropriate handler available")
    hnd$map(subent, extout, rowidcs,
            masterHandler = masterHandler)
  }


  getJac <- function(subent, extout, rowidcs=NULL) {

    hnd <- getResponsibleHandler(subent)
    if (is.null(hnd))
      stop("no appropriate handler available")
    hnd$getJac(subent, extout, rowidcs,
               masterHandler = masterHandler)
  }


  extractData <- function(subent, rowidcs=NULL, ret.values=TRUE) {

    hnd <- getResponsibleHandler(subent)
    if (is.null(hnd))
      stop("no appropriate handler available")
    hnd$extractData(subent, rowidcs, ret.values,
                    masterHandler = masterHandler)
  }


  isHomogeneousMap <- function(subent) {
    hnd <- getResponsibleHandler(subent)
    if (is.null(hnd))
      stop("no appropriate handler available")
    hnd$isHomogeneousMap(subent, masterHandler = masterHandler)
  }


  isLinearMap <- function(subent) {
    hnd <- getResponsibleHandler(subent)
    if (is.null(hnd))
      stop("no appropriate handler available")
    hnd$isLinearMap(subent, masterHandler = masterHandler)
  }


  listHandlers <- function(configured=NA) {

    pickFlag <- if (is.na(configured)) rep(TRUE, length(subentHandlerList))
    else {
      isConfigVec <- sapply(subentHandlerList, function(x) x$isConfigured())
      !xor(isConfigVec, configured)
    }
    sapply(subentHandlerList[pickFlag], function(x) x$getName())
  }


  getHandlerByName <- function(name) {

    for (curHandler in subentHandlerList)
      if (curHandler$getName() == name) {
        return(curHandler)
      }
    NULL
  }


  removeHandlerByName <- function(name) {

    handlerNames <- listHandlers()
    subentHandlerList <<- subentHandlerList[handlerNames != name]
  }


  getResponsibleHandler <- function(subent) {

    for (curHandler in subentHandlerList) {
      if (curHandler$canMap(subent, quiet=TRUE,
                            masterHandler = masterHandler)) {
        return(curHandler)
      }
    }
    NULL
  }


  masterHandler <- list(getName = getName,
                        getInfo = getInfo,
                        # data mapping/manipulation
                        extractData = extractData,
                        canMap = canMap,
                        needs = needs,
                        map = map,
                        getJac = getJac,
                        # info about mapping properties
                        isHomogeneousMap = isHomogeneousMap,
                        isLinearMap = isLinearMap,
                        # additional functionality
                        listHandlers = listHandlers,
                        getHandlerByName = getHandlerByName,
                        removeHandlerByName = removeHandlerByName,
                        getResponsibleHandler = getResponsibleHandler)

  masterHandler
}



#' Create List of Default Subentry Handler
#'
#' @return List of subentry handlers
#' @export
#'
createDefaultSubentHandlerList <- function() {

  handlerList <- list(
    createSubentHandler_ntot(),
    createSubentHandler_ntot_nat()
  )
  names(handlerList) <- sapply(handlerList, function(x) x$getName())
  handlerList
}
