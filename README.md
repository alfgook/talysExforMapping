### talysExforMapping - R package

This package enables the mapping from predictions of TALYS to the observables
recorded in EXFOR entries. For instance, a cross section associated with a
target in natural composition requires the collection and summation of 
data in files produced by several TALYS calculations. 
This package also helps to determine which calculations are necessary to 
predict the observables recorded in EXFOR entries.

## Requirements

The R packages `data.table` and `Matrix` available on CRAN are prerequisites. 
Further, the custom R package `json ExforUtils` is required whose installation
will also be discussed in the next section.

## Installation

The installation can be done from a terminal by executing the commands:
```{bash}
git clone https://github.com/gschnabel/jsonExforUtils.git
R CMD INSTALL jsonExforUtils

git clone https://github.com/gschnabel/talysExforMapping.git
R CMD INSTALL talysExforMapping
```

## Basic usage

First, we load the package and set up some handlers that know how to deal with
EXFOR subentries. 
```{r}
library(talysExforMapping)

subentHandlerList <- createDefaultSubentHandlerList()
subentHandler <- createSubentHandler(subentHandlerList)
exforHandler <- createExforHandler(subentHandler)
```

Next, we load an exemplary EXFOR subentry provided with the package: 
```{r}
xs_ntot_Al27 <- example_subent_10004002
```

If we want to know which TALYS calculations are required to obtain the
prediction corresponding to this EXFOR subentry, we can run
```{r}
expDt <- exforHandler$extractData(list(xs_ntot_Al27), ret.values=FALSE)
exforNeedsDt <- exforHandler$needs(expDt, list(xs_ntot_Al27))
print(exforNeedsDt)
```
The datatable `exforNeedsDt` contains columns indicating for which 
projectile, element, mass combinations TALYS needs to be run for 
a comparison with the EXFOR subentry.
