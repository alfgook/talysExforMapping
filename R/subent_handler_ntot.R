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


#' Create Specific Subent Handler
#'
#' Handles an EXFOR entry with (N,TOT) for a specific target
#'
#' @return Returns closure with functions as given in [.].
#'         In addition, the closure contains the functions \code{setAbuAgent}.
#' @export
#'
createSubentHandler_ntot <- function() {


  getName <- function() {
    "handler_ntot"
  }


  getInfo <- function() {
    ""
  }


  isConfigured <- function() { TRUE }


  configure <- function(config) { TRUE }


  canMap <- function(subent, quiet=TRUE, ret.score=FALSE, masterHandler=NULL) {

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

    # only handling SIG without any extra specifications
    if (!grepl("^,SIG$", reacStruc$quantspec))
      return(retmsg("quantity specification must only contain SIG", 0))

    # not handling natural composition of isotopes
    if (reacStruc$target$A == 0)
      return(retmsg("natural isotope composition not supported", 2))

    # not handling targets in metastable state
    if (!is.null(reacStruc$target$meta))
      return(retmsg("cannot handle target in metastable state", 3))

    # do not accept characters (e.g. MASS) as residual
    if (!is.null(reacStruc$residual) && !is.list(reacStruc$residual))
      return(retmsg("residual must be missing or parsable as nucleus", 3))

    # not handling residual nucleus in metastable state
    if (!is.null(reacStruc$residual$meta))
      return(retmsg("cannot handle residual in metastable state", 3))

    # handling only specific projectile types
    if (!reacStruc$projectile %in% c("N","P","D", "T", "HE3", "A"))
      return(retmsg("can only handle projectile types N, P, D, T, HE3, A", 4))

    particleStr <- generateTalysParticleStr(reacStruc$process)
    # only handling TOT, EL and particle strings
    if (!reacStruc$process %in% c("TOT", "EL", "INL") &&
        is.null(particleStr))
      return(retmsg(paste0("handler expects TOT/EL reaction ",
                           "or a particle string as process"), 4))

    # if particle string, then residual must be specified
    if (!is.null(particleStr) && !is.list(reacStruc$residual))
      return(retmsg(paste0("if process field contains particle string, ",
                           "residual nucleus must be given"), 5))

    # not handling if DATA and EN not present
    if (! all(c("EN", "DATA") %in% subent$DATA$DESCR))
      return(retmsg("EN and DATA column are missing", 6))

    # not handling units different from MEV and MB
    if (subent$DATA$UNIT[subent$DATA$DESCR=="EN"] != "MEV" ||
        subent$DATA$UNIT[subent$DATA$DESCR=="DATA"] != "MB")
      return(retmsg("EN must be in MEV and DATA in MB", 8))

    # missing values in column are not accepted
    if (any(is.na(subent$DATA$TABLE$EN)))
      return(retmsg("missing values in EN are not accepted", 8))

    retmsg("success", 10)
  }


  extractData <- function(subent, rowidcs=NULL, ret.values=TRUE, masterHandler=NULL) {

    stopifnot(canMap(subent, quiet=TRUE))
    if (is.null(rowidcs))
      rowidcs <- seq_len(nrow(subent$DATA$TABLE))

    parsedReac <- parseReacExpr(subent$BIB$REACTION)
    selectedEn <- subent$DATA$TABLE$EN[rowidcs]
    selectedData <- subent$DATA$TABLE$DATA[rowidcs]
    resDt <- data.table(EXPID = subent$ID, DIDX = rowidcs,
                        REAC = parsedReac$reac, L1 = selectedEn, L2 = 0, L3 = 0)
    if (isTRUE(ret.values)) resDt[, DATA := selectedData]
    resDt
  }


  needs <- function(subent, rowidcs=NULL, masterHandler=NULL) {

    if (!canMap(subent, quiet=TRUE))
      stop(paste0("handler ", getName(), " cannot map subent with ID ", subent$ID))
    if (is.null(rowidcs)) rowidcs <- TRUE

    reacStruc <- parseReacExpr(subent$BIB$REACTION)

    proj <- reacStruc$projectile
    targetSym <- reacStruc$target$sym
    targetMass <- reacStruc$target$A
    process <- reacStruc$process
    particleStr <- generateTalysParticleStr(sub("INL", proj, process))  # replace INL by projectile
    energy <- subent$DATA$TABLE$EN[rowidcs]

    talysReacStr <- if (is.null(particleStr)) paste0("CS/", process) else
      paste0("CS/REAC/", particleStr, "/TOT")

    proj <- sub("HE3", "H", proj)  # convert HE3 to talys notation

    unique(data.table(PROJECTILE = proj,
               ELEMENT = targetSym,
               MASS = targetMass,
               REAC = talysReacStr,
               L1 = energy,
               L2 = 0, L3 = 0))
  }


  map <- function(subent, extout, rowidcs=NULL, masterHandler=NULL) {

    if (!canMap(subent))
      stop(paste0("handler ", getName(), " cannot map subent with ID ", subent$ID))
    if (is.null(extout$V1))
      stop("extout must contain a column V1")

    setkeyv(extout, c("PROJECTILE","ELEMENT","MASS","REAC","L1","L2","L3"))

    if (is.null(rowidcs))
      rowidcs <- seq_len(nrow(subent$DATA$TABLE))

    reacStruc <- parseReacExpr(subent$BIB$REACTION)
    curProj <- reacStruc$projectile
    curElem <- reacStruc$target$sym
    curMass <- reacStruc$target$A
    curProc <- reacStruc$process

    particleStr <- generateTalysParticleStr(sub("INL", curProj, curProc))
    talysReacStr <- if (is.null(particleStr)) paste0("CS/", curProc) else
      paste0("CS/REAC/", particleStr, "/TOT")

    curProj <- sub("HE3", "H", curProj)  # convert HE3 to talys notation

    extout[J(curProj,curElem, curMass, talysReacStr), {

      selectedExpEn <- subent$DATA$TABLE$EN[rowidcs]
      modGridRange <- range(L1)
      expGridRange <- range(selectedExpEn)
      if (expGridRange[1] < modGridRange[1] ||
          expGridRange[2] > modGridRange[2])
        stop(paste0("experimental energies outside model grid for subent ID ", subent$ID))

      if (length(L1)==1) {
        if (!all(selectedExpEn == L1))
          stop(paste0("only one energy point in model grid and it does not match experimental energies ",
                      "for subent ID ", subent$ID))
        else
          interpolResult <- as.numeric(V1)
      }
      else  # normal interpolation on grid if length(L1) > 1
        interpolResult <- approx(L1, V1, xout=selectedExpEn, rule=1)$y

      list(EXPID=subent$ID, DIDX=rowidcs, REAC = reacStruc$reac,
           L1 = selectedExpEn, L2 = 0, L3 = 0, DATA = interpolResult)
    }]
  }


  getJac <- function(subent, extout, rowidcs=NULL, masterHandler=NULL) {

    if (!canMap(subent))
      stop(paste0("handler ", getName(), " cannot map subent with ID ", subent$ID))
    if (is.null(extout$IDX))
      stop("extout must contain a column IDX")

    reacStruc <- parseReacExpr(subent$BIB$REACTION)
    curProj <- reacStruc$projectile
    curElem <- reacStruc$target$sym
    curMass <- reacStruc$target$A
    curProc <- reacStruc$process

    particleStr <- generateTalysParticleStr(sub("INL", curProj, curProc))
    talysReacStr <- if (is.null(particleStr)) paste0("CS/", curProc) else
      paste0("CS/REAC/", particleStr, "/TOT")

    curProj <- sub("HE3", "H", curProj)  # convert HE3 to talys notation

    if (is.null(rowidcs))
      rowidcs <- seq_len(nrow(subent$DATA$TABLE))

    setkeyv(extout, c("PROJECTILE","ELEMENT","MASS","REAC","L1","L2","L3"))
    modRes <- extout[J(curProj,curElem,curMass, talysReacStr), list(IDX, L1)]

    selectedExpEn <- subent$DATA$TABLE$EN[rowidcs]
    modGridRange <- range(modRes$L1)
    expGridRange <- range(selectedExpEn)
    if (expGridRange[1] < modGridRange[1] ||
        expGridRange[2] > modGridRange[2])
      stop(paste0("experimental energies outside model grid for subent ID ", subent$ID))

    if (length(modRes$L1)==1) {
      if (!all(selectedExpEn == modRes$L1))
        stop(paste0("only one energy point in model grid and it does not match experimental energies ",
                    "for subent ID ", subent$ID))
      else
        sparseMatrix(i=1, j=modRes$IDX, x=1, dims=c(1, nrow(extout)))
    }
    else  {  # normal interpolation on grid if length(L1) > 1

      idx1 <- findInterval(selectedExpEn, modRes$L1)  # L1 is sorted!

      stopifnot(idx1 > 0, idx1 <= length(modRes$L1))
      isOnEdge <- idx1==length(modRes$L1)
      idx1[isOnEdge] <- idx1[isOnEdge] - 1  # special case if point on right edge

      idx2 <- idx1 + 1
      en1 <- modRes$L1[idx1]
      en2 <- modRes$L1[idx2]
      len <- en2 - en1

      rowi <- rep(seq_along(rowidcs), 2)
      colj <- c(idx1, idx2)  # colj is only with respect to the selected rows
      # has to be converted to the global index in extout later
      valx <- c((en2 - selectedExpEn) / len,
                (selectedExpEn - en1) / len)
      rowi <- rowi[valx > .Machine$double.eps * 100]
      colj <- colj[valx > .Machine$double.eps * 100]
      valx <- valx[valx > .Machine$double.eps * 100]

      sparseMatrix(i=rowi, j=modRes$IDX[colj], x=valx,
                   dims=c(length(rowidcs), nrow(extout)))
    }
  }


  isLinearMap <- function(subent, masterHandler = NULL) { TRUE }

  isHomogeneousMap <- function(subent, masterHandler = NULL) { TRUE }


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
