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

#' Create Specific Subentry Handler
#'
#' Handles an EXFOR entry with (N,TOT) and target in natural composition
#'
#' @return Returns closure with functions as given in [.].
#'         In addition, the closure contains the functions \code{setAbuAgent}.
#' @export
#'
createSubentHandler_ntot_nat = function() {

  getName <- function() {
    "handler_ntot_nat"
  }


  getInfo <- function() {
    ""
  }


  isConfigured <- function() {

    !is.null(abuAgent)
  }


  configure <- function(config) {

    if (!is.list(config))
      stop("config must be a list")
    if (is.null(config$abuAgent))
      stop("expecting an abuAgent in config")

    abuAgent <<- config$abuAgent
  }


  canMap <- function(subent, quiet=TRUE, ret.score=FALSE,
                     masterHandler=NULL) {

    retmsg <- function(msg, score) {
      if (!quiet) cat(msg, "\n"); if (ret.score) score else score >= 10
    }

    # not handling different observables in the same table
    if (length(subent$BIB$REACTION) != 1)
      return(retmsg("BIB.REACTION must be a vector of length 1", 0))

    reacStruc <- parseReacExpr(subent$BIB$REACTION, quiet=quiet)
    # check if reaction string can be parsed
    if (is.null(reacStruc))
      return(retmsg("cannot parse BIB.REACTION", 0))

    # not handling ratios, sums of cross sections etc.
    if (reacStruc$type != "reac")
      return(retmsg("arithmetic expression of reactions not supported", 0))

    # only handles natural composition of isotopes
    if (reacStruc$target$A != 0)
      return(retmsg("expects isotope in natural composition", 2))

    # masterHandler must be provided with handler_ntot
    if (is.null(masterHandler))
      return(retmsg("masterHandler must be provided", 9))

    # abuAgent must be created and passed to this handler
    if (is.null(abuAgent))
      return(retmsg("abuAgent must be specified before using handler_ntot_nat"))

    abuDt <-abuAgent$getAbundances(reacStruc$target$sym)
    # handler_ntot must be able to map specific isotopes
    if (is.null(abuDt))
      return(retmsg("abuAgent could not provide data needed on abundances", 9))

    reacStruc$target$A <- abuDt$MASS[1]
    stopifnot(reacStruc$target$A != 0)
    reacStruc$residual <- determineResidualNucleus(reacStruc$projectile,
                                                   reacStruc$target,
                                                   reacStruc$process)

    # residual nucleus could not be determined
    if (is.null(reacStruc$residual) && ! reacStruc$process %in% c("EL", "TOT")) {
      return(retmsg("could not determine residual nucleus", 9))
    }

    reacStr <- reacStrucToStr(reacStruc)
    subent$BIB$REACTION <- reacStr

    if (!masterHandler$canMap(subent, quiet=TRUE)) {
      if (!quiet) {
        cat("### start output of master handler ###\n")
        masterHandler$canMap(subent, quiet)
        cat("### end output of master handler ###\n")
      }
      return(retmsg("no handler available to map specific isotopes", 9))
    }

    retmsg("success", 10)
  }


  needs <- function(subent, rowidcs=NULL, masterHandler=NULL) {

    stopifnot(canMap(subent, quiet=TRUE, masterHandler = masterHandler))

    if (is.null(rowidcs)) rowidcs <- TRUE
    reacStruc <- parseReacExpr(subent$BIB$REACTION)

    targetSym <- reacStruc$target$sym
    abuDt <- abuAgent$getAbundances(targetSym)
    stopifnot(all(abuDt$MASS != 0))

    unique(rbindlist(mapply(function(A) {
      reacStruc$target$A <- A
      reacStruc$residual <- determineResidualNucleus(reacStruc$projectile,
                                                     reacStruc$target,
                                                     reacStruc$process)
      reacStr <- reacStrucToStr(reacStruc)
      subent$BIB$REACTION <- reacStr
      masterHandler$needs(subent, rowidcs)
    }, A=abuDt$MASS, SIMPLIFY = FALSE)))
  }


  extractData <- function(subent, rowidcs=NULL, ret.values=TRUE,
                          masterHandler=NULL) {

    stopifnot(canMap(subent, quiet=TRUE, masterHandler = masterHandler))
    if (is.null(rowidcs))
      rowidcs <- seq_len(nrow(subent$DATA$TABLE))

    parsedReac <- parseReacExpr(subent$BIB$REACTION)
    selectedEn <- subent$DATA$TABLE$EN[rowidcs]
    selectedData <- subent$DATA$TABLE$DATA[rowidcs]
    resDt <- data.table(EXPID = subent$ID, DIDX = rowidcs,
                        REAC = parsedReac$reac, L1 = selectedEn, L2 = 0, L3 = 0)
    if (isTRUE(ret.values)) resDt[, DATA := selectedData]
    resDt[]
  }


  map <- function(subent, extout, rowidcs=NULL, masterHandler=NULL) {

    if (!canMap(subent, masterHandler = masterHandler))
      stop(paste0("handler ", getName(), " cannot map subent with ID ", subent$ID))
    if (is.null(extout$V1))
      stop("extout must contain a column V1")

    setkeyv(extout, c("PROJECTILE","ELEMENT","MASS","REAC","L1","L2","L3"))

    if (is.null(rowidcs))
      rowidcs <- seq_len(nrow(subent$DATA$TABLE))

    reacStruc <- parseReacExpr(subent$BIB$REACTION)
    targetSym <- reacStruc$target$sym
    abuDt <- abuAgent$getAbundances(targetSym)

    # dirty hack? change BIB$REACTION and invoke handler_ntot
    resData <- 0
    for (curRow in seq_len(nrow(abuDt))) {
      reacStruc$target$A <- abuDt[curRow, MASS]
      reacStruc$residual <- determineResidualNucleus(reacStruc$projectile,
                                                     reacStruc$target,
                                                     reacStruc$process)
      subent$BIB$REACTION <- reacStrucToStr(reacStruc)
      curWeight <- abuDt[curRow, ABU]
      curMap <- masterHandler$map(subent, extout, rowidcs)
      curData <- curMap$DATA * curWeight
      resData <- resData + curData
    }
    resDt <- curMap[, REAC:=reacStruc$reac]
    resDt <- curMap[, DATA:=resData]
    resDt[]
  }


  getJac <- function(subent, extout, rowidcs=NULL, masterHandler=NULL) {

    if (!canMap(subent, masterHandler = masterHandler))
      stop(paste0("handler ", getName(), " cannot map subent with ID ", subent$ID))
    if (is.null(extout$V1))
      stop("extout must contain a column V1")

    setkeyv(extout, c("PROJECTILE","ELEMENT","MASS","REAC","L1","L2","L3"))

    if (is.null(rowidcs))
      rowidcs <- seq_len(nrow(subent$DATA$TABLE))

    reacStruc <- parseReacExpr(subent$BIB$REACTION)
    targetSym <- reacStruc$target$sym
    abuDt <- abuAgent$getAbundances(targetSym)

    # dirty hack? change BIB$REACTION and invoke handler_ntot
    resData <- 0
    resDf <- NULL
    for (curRow in seq_len(nrow(abuDt))) {
      reacStruc$target$A <- abuDt[curRow, MASS]
      reacStruc$residual <- determineResidualNucleus(reacStruc$projectile,
                                                     reacStruc$target,
                                                     reacStruc$process)
      subent$BIB$REACTION <- reacStrucToStr(reacStruc)
      curWeight <- abuDt[curRow, ABU]
      curJac <- masterHandler$getJac(subent, extout, rowidcs)
      curDf <- summary(curJac)
      curDf$x <- curDf$x * curWeight
      resDf <- rbind(resDf, curDf)
    }
    sparseMatrix(i=resDf$i, j=resDf$j, x=resDf$x, dims=c(length(rowidcs), nrow(extout)))
  }


  isLinearMap <- function(subent, masterHandler = NULL) {

    if (!canMap(subent, masterHandler = masterHandler))
      stop(paste0("handler ", getName(), " cannot map subent with ID ", subent$ID))

    reacStruc <- parseReacExpr(subent$BIB$REACTION)
    targetSym <- reacStruc$target$sym
    abuDt <- abuAgent$getAbundances(targetSym)
    isLinearDt <- abuDt[, list(isLinear = {
      reacStruc$target$A <- MASS
      reacStruc$residual <- determineResidualNucleus(reacStruc$projectile,
                                                     reacStruc$target,
                                                     reacStruc$process)
      subent$BIB$REACTION <- reacStrucToStr(reacStruc)
      masterHandler$isLinearMap(subent)
    }), by=c("CHARGE", "MASS")]
    all(isLinearDt[, isLinear])
  }


  isHomogeneousMap <- function(subent, masterHandler = NULL) {

    if (!canMap(subent, masterHandler = masterHandler))
      stop(paste0("handler ", getName(), " cannot map subent with ID ", subent$ID))

    reacStruc <- parseReacExpr(subent$BIB$REACTION)
    targetSym <- reacStruc$target$sym
    abuDt <- abuAgent$getAbundances(targetSym)
    isHomogeneousDt <- abuDt[, list(isHomogeneous = {
      reacStruc$target$A <- MASS
      reacStruc$residual <- determineResidualNucleus(reacStruc$projectile,
                                                     reacStruc$target,
                                                     reacStruc$process)
      subent$BIB$REACTION <- reacStrucToStr(reacStruc)
      masterHandler$isHomogeneousMap(subent)
    }), by=c("CHARGE", "MASS")]
    all(isHomogeneousDt[, isHomogeneous])
  }

  # contains config

  abuAgent <- NULL


  list(getName = getName,
       getInfo = getInfo,
       # special configuation
       configure = configure,
       isConfigured = isConfigured,
       # data mapping/manipulation
       canMap = canMap,
       extractData = extractData,
       needs = needs,
       map = map,
       getJac = getJac,
       # mapping properties
       isLinearMap = isLinearMap,
       isHomogeneousMap = isHomogeneousMap)
}



