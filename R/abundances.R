

#' Create Abundance Agent
#'
#' Constructs a closure with function to set and query abundances of elements.
#'
#' @param talysStructurePath path to the structure/abundance directory of TALYS
#'
#' @return List with functions
#'         \tabular{ll}{
#'           \code{setAbundances} \tab TODO \cr
#'           \code{getAbundances(sym)} \tab
#'           get a data table with abundances of isotopes of an element \code{sym}, e.g. 'Fe'
#'         }
#'
#' @export
#'
createAbuAgent <- function(talysStructurePath) {


  normalizeSym <- function(sym) {
    stopifnot(nchar(sym)==1 || nchar(sym)==2)
    paste0(toupper(substring(sym,1,1)),
           tolower(substring(sym,2,2)))
  }

  setAbundances <- function(sym) {
    stop("TODO")
  }


  getAbundances <- function(sym) {

    stopifnot(length(sym)==1)
    sym <- normalizeSym(sym)
    filename <- paste0(sym, ".abun")
    filepath <- file.path(talysStructurePath, filename)
    if (!file.exists(filepath)) return(NULL)

    abuDt <- try(read.table(filepath, colClasses=c(rep("integer", 2),
                                               rep("numeric", 2),
                                               "character"),
                        col.names = c("CHARGE", "MASS", "ABU", "XXX", "ELEMENT"),
                        stringsAsFactors = FALSE))
    if ("try-error" %in% class(abuDt)) return(NULL)

    abuDt <- as.data.table(abuDt)
    abuDt$ELEMENT <- sub("^[0-9]+","", abuDt$ELEMENT)
    abuDt$ABU <- abuDt$ABU / 100
    abuDt
  }

  list(getAbundances=getAbundances)
}
