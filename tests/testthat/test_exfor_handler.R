context("Test EXFOR Handler")
library(talysExforMapping)
library(MongoEXFOR)

subentHandlerList <- createDefaultSubentHandlerList()
subentHandler <- createSubentHandler(subentHandlerList)
exforHandler <- createExforHandler(subentHandler)
# subentList <- replicate(100, example_subent_10004002, simplify = FALSE)
subentList <- c(example_subents_ntot, example_subents_ntot_nat)

abuAgent <- createAbuAgent("/home/georg/bin/talys/structure/abundance/")
subentHandler$getHandlerByName("handler_ntot_nat")$configure(list(abuAgent=abuAgent))


test_that("exforHandler determines correctly mapable subentries",{

  expect_true(all(exforHandler$canMap(subentList)))
})


test_that("exforHandler$needs returns data table with correct structure", {

  expDt <- exforHandler$extractData(subentList[1], ret.values=FALSE)
  exforNeedsDt <- exforHandler$needs(expDt, subentList[1])
  subentNeedsDt <- subentHandler$needs(subentList[[1]])
  expect_true("IDX" %in% names(exforNeedsDt))
  exforNeedsDt[, IDX:=NULL]
  setkey(exforNeedsDt, NULL)
  expect_identical(exforNeedsDt, subentNeedsDt)
})


test_that("exforHandler$extractData returns data table with correct structure", {

  exforExtout <- exforHandler$extractData(subentList[1])
  subentExtout <- subentHandler$extractData(subentList[[1]])
  expect_true("IDX" %in% names(exforExtout))
  exforExtout[, IDX:=NULL]
  expect_identical(exforExtout, subentExtout)
})


test_that("exforHandler$needs is consistent with subentHandler$needs", {

  for (i in seq_along(subentList)) {
    curSubent <- subentList[[i]]
    expDt <- exforHandler$extractData(list(curSubent), ret.values=FALSE)
    subentNeedsDt <- subentHandler$needs(curSubent)
    exforNeedsDt <- exforHandler$needs(expDt, list(curSubent))
    exforNeedsDt[,IDX:=NULL]
    setkey(exforNeedsDt, PROJECTILE, ELEMENT, MASS, REAC, L1, L2, L3)
    setkey(subentNeedsDt, PROJECTILE, ELEMENT, MASS, REAC, L1, L2, L3)
    expect_identical(exforNeedsDt, subentNeedsDt)
  }
})

test_that("exforHandler$extractData consistent with subentHandler$extractData", {

  for (i in seq_along(subentList)) {

    curSubent <- subentList[[i]]
    expDt <- exforHandler$extractData(list(curSubent), ret.values=FALSE)
    exforNeedsDt <- exforHandler$needs(expDt, list(curSubent))
    exforNeedsDt[, V1:=seq_len(.N)]
    exforExpDt <- exforHandler$extractData(list(curSubent))

    subentExpDt <- subentHandler$extractData(curSubent)
    subentMappedDt <- subentHandler$map(curSubent, exforNeedsDt)
    exforMappedDt <- exforHandler$map(exforExpDt, exforNeedsDt, list(curSubent))
    setkey(exforMappedDt, EXPID, DIDX, REAC, L1, L2, L3)
    setkey(subentMappedDt, EXPID, DIDX, REAC, L1, L2, L3)
    exforMappedDt[, IDX:=NULL]

    expect_identical(exforMappedDt, subentMappedDt,
                     info = paste0("problem with subent ", curSubent$ID,
                                   " - loop counter is ", i))
  }
})


test_that("exforHandler$extractData consistent with subentHandler$extractData when using rowidcs", {

  expDt <- exforHandler$extractData(subentList[1], ret.values=FALSE)
  exforNeedsDt <- exforHandler$needs(expDt, subentList[1])
  exforNeedsDt[, V1:=seq_len(.N)]
  exforExpDt <- exforHandler$extractData(subentList[1])
  exforExpDt <- exforExpDt[DIDX >=10 & DIDX <= 20,]

  subentNeedsDt <- subentHandler$needs(subentList[[1]])
  subentNeedsDt[, V1:=seq_len(.N)]

  perm <- sample(10:20, 11)
  exforExpDt[, DIDX:=perm]
  subentMappedDt <- subentHandler$map(subentList[[1]], subentNeedsDt, perm)
  exforMappedDt <- exforHandler$map(exforExpDt, exforNeedsDt, subentList[1])

  setkey(exforMappedDt, IDX)
  exforMappedDt[, IDX:=NULL]
  setkey(exforMappedDt, NULL)
  expect_identical(exforMappedDt, subentMappedDt)
})


test_that("exforHandler$map yields the same result as application of exforHandler$getJac", {

  expDt <- exforHandler$extractData(subentList, ret.values=FALSE)
  exforNeedsDt <- exforHandler$needs(expDt, subentList)
  exforNeedsDt[, V1:=seq_len(.N)]

  mapRes1 <- exforHandler$map(expDt, exforNeedsDt, subentList)[,DATA[order(IDX)]]
  jac <- exforHandler$getJac(expDt, exforNeedsDt, subentList)
  setkey(expDt, IDX)
  setkey(exforNeedsDt, IDX)
  mapRes2 <- as.numeric(jac %*% exforNeedsDt[,V1])
  expect_identical(mapRes1, mapRes2)
})


test_that("exforHandler$map yields same result as exforHandler$getJac with filtered and permuted rows", {

  expDt <- exforHandler$extractData(subentList, ret.values=FALSE)
  exforNeedsDt <- exforHandler$needs(expDt, subentList)
  exforNeedsDt[, V1:=seq_len(.N)]

  filteredPerm <- sample.int(nrow(expDt), 30)
  expDt <- expDt[filteredPerm,]
  expDt[, IDX:=seq_len(.N)]
  mapRes1 <- exforHandler$map(expDt, exforNeedsDt, subentList)
  setkey(mapRes1, IDX)
  mapRes1 <- mapRes1[,DATA]
  setkey(expDt, IDX)
  setkey(exforNeedsDt, IDX)
  jac <- exforHandler$getJac(expDt, exforNeedsDt, subentList)
  mapRes2 <- as.numeric(jac %*% exforNeedsDt[,V1])
  expect_identical(mapRes1, mapRes2)
})


test_that("exforHandler$map provides correct result if duplicates in expDt", {

  expDt <- exforHandler$extractData(subentList[1], ret.values=FALSE)
  expDt <- expDt[1:5,]
  expDt2 <- rbind(expDt, expDt)
  expDt2[, IDX:=seq_len(.N)]
  exforNeedsDt <- exforHandler$needs(expDt, subentList[1])
  exforNeedsDt[, V1:=seq_len(.N)]
  mapRes1 <- exforHandler$map(expDt, exforNeedsDt, subentList[1])
  mapRes2 <- exforHandler$map(expDt2, exforNeedsDt, subentList[1])
  setkey(mapRes1, IDX)
  doubleMapRes1 <- rbind(mapRes1, mapRes1)
  doubleMapRes1[, IDX:=seq_len(.N)]
  setkey(doubleMapRes1, IDX)
  setkey(mapRes2, IDX)
  expect_identical(doubleMapRes1, mapRes2)
})

# subentHandlerList <- createDefaultSubentHandlerList()
# subentHandler <- createSubentHandler(subentHandlerList)
# exforHandler <- createExforHandler(subentHandler)
# subentList <- replicate(200, example_subent_10004002, simplify = FALSE)
#
#
# expDt <- exforHandler$extractData(subentList, ret.values=FALSE)
# system.time(needsDt <- exforHandler$needs(expDt, subentList))
#
# needsDt[,V1:=seq_len(.N)]
#
# system.time(mapRes <- exforHandler$map(expDt, needsDt, subentList))
# system.time(sensMat <- exforHandler$getJac(expDt, needsDt, subentList))
#
# sensMat2 <- subentHandler$getJac(subentList[[1]], needsDt, rowidcs = 5:3)
#
# (sensMat %*% needsDt[,V1])[1:10,]
# (t(sensMat) %*% expDt[,DATA])[1:10,]
