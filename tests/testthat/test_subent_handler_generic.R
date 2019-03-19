context("generic tests for all subent handlers")
library(talysExforMapping)
library(jsonExforUtils)
library(data.table)

subentHandlerList <- createDefaultSubentHandlerList()
subentHandler <- createSubentHandler(subentHandlerList)
exforHandler <- createExforHandler(subentHandler)

# initialize the abundance agent
abuAgent <- createAbuAgent("/home/georg/bin/talys/structure/abundance/")
subentHandler$getHandlerByName("handler_ntot_nat")$configure(list(abuAgent = abuAgent))
#subentHandler$removeHandlerByName("handler_ntot")


test_that("individual handlers decide correctly what they can map", {

  masterHandler <- subentHandler
  handlerNames <- subentHandler$listHandlers()
  for (curName in handlerNames) {
    hnd <- subentHandler$getHandlerByName(curName)
    curSubents <- get(sub("handler_", "example_subents_", curName))
    canMapFlag <- sapply(curSubents, hnd$canMap, masterHandler=masterHandler, quiet=FALSE)
    expect_true(all(canMapFlag),
                info=paste0("handler ", hnd$getName(), " cannot map subentries with ID ",
                             paste0(sapply(curSubents[!canMapFlag], function(x) x$ID), collapse=", "),
                            " and reactions ", paste0(sapply(curSubents[!canMapFlag], function(x) x$BIB$REACTION), collapse=", ")))
  }
})


test_that("method 'needs' of all handlers yields a data table with correct structure", {

  masterHandler <- subentHandler
  handlerNames <- subentHandler$listHandlers()
  for (curName in handlerNames) {
    hnd <- subentHandler$getHandlerByName(curName)
    curSubents <- get(sub("handler_", "example_subents_", curName))
    needsDt <- hnd$needs(curSubents[[1]], masterHandler = subentHandler)
    expect_true(all(c("PROJECTILE", "ELEMENT", "MASS", "REAC", "L1", "L2", "L3") %in%
                    names(needsDt)),
                info = paste0("method 'needs' of handler ", hnd$getName(),
                               " does not yield a correct data table for subentry with ID",
                               curSubents[[1]]$ID))
  }
})


test_that("method 'extractData' of all handlers returns correct structure", {

  masterHandler <- subentHandler
  handlerNames <- subentHandler$listHandlers()
  for (curName in handlerNames) {
    hnd <- subentHandler$getHandlerByName(curName)
    curSubents <- get(sub("handler_", "example_subents_", curName))
    extout <- hnd$extractData(curSubents[[1]], masterHandler = subentHandler)
    expect_true(all(c("EXPID", "DIDX", "DATA") %in% names(extout)),
                info = paste0("method 'extractData' of handler ", hnd$getName(),
                               " does not yield a correct data table for subentry with ID",
                               curSubents[[1]]$ID))
  }
})


test_that("method 'extractData' of all handlers returns correct number of rows", {

  masterHandler <- subentHandler
  handlerNames <- subentHandler$listHandlers()
  for (curName in handlerNames) {
    hnd <- subentHandler$getHandlerByName(curName)
    curSubents <- get(sub("handler_", "example_subents_", curName))
    extout <- hnd$extractData(curSubents[[1]], rowidcs=2:5, masterHandler = subentHandler)
    expect_equal(nrow(extout), 4,
                info = paste0("method 'extractData' of handler ", hnd$getName(),
                               " does not yield a data table with correct number of rows ",
                              "for subentry with ID", curSubents[[1]]$ID))
  }
})


test_that("method 'map' and 'getJac' of all handlers are consistent with each other", {

  masterHandler <- subentHandler
  handlerNames <- subentHandler$listHandlers()
  for (curName in handlerNames) {
    hnd <- subentHandler$getHandlerByName(curName)
    curSubents <- get(sub("handler_", "example_subents_", curName))
    if (! hnd$isHomogeneousMap(curSubents[[1]], masterHandler = masterHandler)) next
    needsDt <- hnd$needs(curSubents[[1]], masterHandler = masterHandler)
    needsDt[, V1 := seq_len(.N)]
    needsDt[, IDX := seq_len(.N)]
    selIdcs <- sample.int(10,5)
    mapRes1 <- hnd$map(curSubents[[1]], needsDt, selIdcs, masterHandler = masterHandler)[,DATA]
    jac <- hnd$getJac(curSubents[[1]], needsDt, selIdcs, masterHandler = masterHandler)
    setkey(needsDt, IDX)
    mapRes2 <- as.numeric(jac %*% needsDt[,V1])
    expect_equal(mapRes1, mapRes2,
                 info = paste0("method 'map' and 'getJac' of handler ", hnd$getName(),
                               "are inconsistent with each other ",
                               "for subentry with ID", curSubents[[1]]$ID))
  }
})






