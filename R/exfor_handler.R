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

#' Create EXFOR Handler
#'
#' @return Returns an EXFOR handler
#' @export
#' @import data.table Matrix jsonExforUtils
#'
createExforHandler <- function(subentHandler) {

  # structure specification of data tables

  reqColsInNeedsDt <- c("IDX", "PROJECTILE", "ELEMENT", "MASS",
                        "REAC", "L1", "L2","L3")

  reqColsInExtout <- c(reqColsInNeedsDt, "V1")

  reqColsInExpdt <- c("IDX", "EXPID", "DIDX")

  # public functions

  canMap <- function(subents, quiet=TRUE) {

    if (!is.list(subents) || !all(sapply(subents, is.list)))
      stop(paste0("subents must be a list of subents"))

    mapable <- sapply(subents, subentHandler$canMap, quiet=TRUE)
    if (!all(mapable) && !quiet) {
      mapply(function(i,x) {
        cat("=== IDX ", i, " ===\n")
        subentHandler$canMap(x, quiet=FALSE)
      }, i=which(!mapable), x=subents[!mapable])
    }
    mapable
  }

  needs <- function(expDt, subents) {

    if (!is.list(subents) || !all(sapply(subents, is.list)))
      stop(paste0("subents must be a list of subents"))
    if (!all(reqColsInExpdt %in% names(expDt)))
      stop(paste0("expdt must contain columns ",
                  paste0(reqColsInExpdt, collapse=", ")))

    setkey(expDt, EXPID)

    subentDt <- data.table(EXPID = sapply(subents, function(x) x$ID),
                               LISTPOS = seq_along(subents), key="EXPID")
    if (anyDuplicated(subentDt$EXPID))
      stop("No duplicates within the subents list are permitted")

    joinedDt <- merge(expDt, subentDt, all.x=TRUE)
    if (any(is.na(joinedDt[,LISTPOS])))
      stop(paste0("EXPIDs ", paste0(head(joinedDt[is.na(LISTPOS), unique(EXPID)], 3), collapse=", "),
                  " are not present in subents"))

    joinedDt <- joinedDt[, subentHandler$needs(subents[[LISTPOS[1]]], DIDX),by="EXPID"]

    keepNames <- c("PROJECTILE", "ELEMENT", "MASS", "REAC", "L1", "L2", "L3")
    dismissNames <- setdiff(names(joinedDt), keepNames)
    joinedDt[, (dismissNames) := NULL]
    setkeyv(joinedDt, keepNames)
    joinedDt <- unique(joinedDt, by=keepNames)
    joinedDt[, IDX:=seq_len(.N)]
    setcolorder(joinedDt, c("IDX", keepNames))
    joinedDt[]
  }

  map <- function(expDt, needsDt, subents) {

    if (!is.list(subents) || !all(sapply(subents, is.list)))
      stop(paste0("subents must be a list of subents"))
    if (!all(reqColsInExtout %in% names(needsDt)))
      stop(paste0("extout must contain columns ",
                  paste0(reqColsInExtout, collapse=", ")))

    if (!all(reqColsInExpdt %in% names(expDt)))
      stop(paste0("expdt must contain columns ",
                  paste0(reqColsInExpdt, collapse=", ")))
    if (any(expDt[,is.na(EXPID) | is.na(DIDX)]))
      stop("expdt is not allowed to have NA values in columns EXPID, DIDX")

    subentDt <- data.table(EXPID = sapply(subents, function(x) x$ID),
                           LISTPOS = seq_along(subents), key="EXPID")
    subentDt <- unique(subentDt, by="EXPID")

    setkey(expDt, EXPID, DIDX)
    joinedDt <- merge(expDt, subentDt)

    helperFun <- function(x, y, z) {
      subentHandler$map(x, y, z)[, list(DIDX, REAC, L1, L2, L3, DATA)]
    }

    joinedDt[, c("DIDX", "REAC", "L1", "L2", "L3", "DATA") :=
               helperFun(subents[[LISTPOS[1]]], needsDt, DIDX), by="EXPID"]
    joinedDt[, LISTPOS := NULL]
    joinedDt[]
  }


  getJac <- function(expDt, needsDt, subents) {

    if (!is.list(subents) || !all(sapply(subents, is.list)))
      stop(paste0("subents must be a list of subents"))
    if (!all(reqColsInNeedsDt %in% names(needsDt)))
      stop(paste0("needsDt must contain columns ",
                  paste0(reqColsInExtout, collapse=", ")))
    if (!all(reqColsInExpdt %in% names(expDt)))
      stop(paste0("expdt must contain columns ",
                  paste0(reqColsInExpdt, collapse=", ")))
    if (any(expDt[,is.na(EXPID) | is.na(DIDX)]))
      stop("expdt is not allowed to have NA values in columns EXPID, DIDX")
    if (min(expDt$IDX) < 1 || max(expDt$IDX > nrow(expDt)))
      stop("some IDX in expDt is out of range")
    if (min(needsDt$IDX) < 1 || max(needsDt$IDX) > nrow(needsDt))
      stop("some IDX in needsDt is out of range")

    setkey(expDt, EXPID)
    subentDt <- data.table(EXPID = sapply(subents, function(x) x$ID),
                           LISTPOS = seq_along(subents), key="EXPID")
    subentDt <- unique(subentDt, by="EXPID")
    joinedDt <- merge(expDt, subentDt)

    spec <- joinedDt[, {
      subMat <- subentHandler$getJac(subents[[LISTPOS[1]]], needsDt, DIDX)
      df <- summary(subMat)
      df$i <- IDX[df$i]
      df[df$x > .Machine$double.eps*100, ]
    } ,by="EXPID"]

    sparseMatrix(i=spec$i, j=spec$j, x=spec$x,
                                   dims=c(nrow(expDt), nrow(needsDt)))
  }


  extractData <- function(subents, expDt=NULL, ret.values=TRUE) {

    if (!is.list(subents) || !all(sapply(subents, is.list)))
      stop(paste0("subents must be a list of subents"))
    if (is.null(expDt)) {
      resDt <- rbindlist(lapply(subents, subentHandler$extractData,
                             ret.values=ret.values))
      oldNames <- copy(names(resDt))
      resDt[,IDX:=seq_len(.N)]
      setcolorder(resDt, c("IDX", oldNames))
    }
    else
    {
      if (!all(reqColsInExpdt %in% names(expDt)))
        stop(paste0("expdt must contain columns ",
                    paste0(reqColsInExpdt, collapse=", ")))
      if (any(expDt[,is.na(EXPID) | is.na(DIDX)]))
        stop("expdt is not allowed to have NA values in columns EXPID, DIDX")

      setkey(expDt, EXPID, DIDX)

      resDt <- expDt
      if (ret.values)
        for (subent in subents) {
          expDt[J(subent$ID), DATA:=subentHandler$extractData(subent, DIDX, ret.values)$DATA]
        }
    }

    resDt[]
  }


  list(canMap = canMap,
       needs = needs,
       map = map,
       getJac = getJac,
       extractData = extractData)
}





