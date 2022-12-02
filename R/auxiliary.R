
projZmap <- c(N=0, P=1, D=1, T=1, HE3=2, A=2)
projAmap <- c(N=1, P=1, D=2, T=3, HE3=3, A=4)

Zmap <- c(H=1,HE=2,LI=3,BE=4,B=5,C=6,N=7,O=8,F=9,NE=10,
          `NA`=11,MG=12,AL=13,SI=14,P=15,S=16,CL=17,AR=18,  # tricky, NA ambiguous: Natrium or Not Answered?
          K=19,CA=20,SC=21,TI=22,V=23,CR=24,MN=25,FE=26,CO=27,NI=28,CU=29,
          ZN=30,GA=31,GE=32,AS=33,SE=34,BR=35,KR=36,
          RB=37,SR=38,Y=39,ZR=40,NB=41,MO=42,TC=43,RU=44,RH=45,PD=46,AG=47,
          CD=48,IN=49,SN=50,SB=51,TE=52,I=53,XE=54,
          CS=55,BA=56,LA=57,CE=58,PR=59,ND=60,PM=61,SM=62,EU=63,GD=64,
          TB=65,DY=66,HO=67,ER=68,TM=69,YB=70,LU=71,
          HF=72,TA=73,W=74,RE=75,OS=76,IR=77,PT=78,AU=79,HG=80,
          TL=81,PB=82,BI=83,PO=84,AT=85,RN=86,
          FR=87,RA=88,AC=89,TH=90,PA=91,U=92,NP=93,PU=94,AM=95,CM=96,
          BK=97,CF=98,ES=99,FM=100,MD=101,NO=102,LR=103,
          RF=104,DB=105,SG=106,BH=107,HS=108,MT=109,DS=110,RG=111,
          CN=112,NH=113,FL=114,MC=115,LV=116,TS=117,OG=118)


generateTalysParticleStr <- function(processStr) {

  if (length(processStr) != 1)
    stop("processStr must be of length 1")

  processStr <- toupper(processStr)
  processStr <- gsub(" ","", processStr)
  numXpart <- strsplit(processStr, "+", fixed = TRUE)[[1]]
  pat <- "^([0-9]?)([NPDTA]|HE3)$"
  if (all(grepl(pat, numXpart))) {

    regRes <- regexec(pat, numXpart)
    regStrs <- regmatches(numXpart, regRes)
    particleStrs <- sapply(regStrs, function(x) x[3])
    numStrs <- sapply(regStrs, function(x) x[2])
    numAry <- rep(0, 6)
    names(numAry) <- c("N","P","D","T","HE3","A")
    numAry[particleStrs] <- ifelse (numStrs=="", 1, as.numeric(numStrs))
    paste0(numAry, collapse="")
  } else if(processStr=="X") {
    "XXXXXX"
  } else NULL
}


determineResidualNucleus <- function(proj, target, process) {

  tarZ <- target$Z
  tarA <- target$A
  talysParticleStr <- generateTalysParticleStr(process)
  if (! proj %in% c("N", "P", "D", "T", "HE3", "A")) return(NULL)
  if (is.null(talysParticleStr)) return(NULL)

  stopifnot(is.null(target$meta))
  stopifnot(length(tarZ)==1, length(tarA)==1, length(proj)==1)
  pN <- as.integer(strsplit(talysParticleStr, "")[[1]])
  stopifnot(length(pN)==6)
  residualZ <- tarZ - sum(pN * projZmap) + projZmap[proj]
  residualA <- tarA - sum(pN * projAmap) + projAmap[proj]
  stopifnot(!is.na(residualZ), !is.na(residualA))
  sym <- names(Zmap)[Zmap == residualZ[1]]
  stopifnot(length(sym)==1)
  list(Z = residualZ,
       A = residualA,
       sym = sym)
}

