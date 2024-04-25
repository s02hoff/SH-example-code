#######################
#creating capture hist#
#######################
rm(list=ls())

#load libraries
library(plotrix) #need this for the pasteCols command
library(BBmisc) #to standardize covs
library(corrplot) #test correlations of covs
# set wd


### Read in detection data and covariates
LIdet <- read.csv("LI17det_prop.csv", skip = 1, na.strings = 'NA')
options(width=100)
print(head(LIdet))
str(LIdet)
sitecovs <- read.csv("SiteCov_test.csv", skip = 1, na.strings = 'NA')
#sitecovs[is.na(sitecovs)] <- 999 #putting in an obvious incorrect value for NAs because RMark won't allow missing covariate values

# detection history
dethist <- LIdet[,2:25] # 24 visits
dethist[is.na(dethist)] <- "." # in RMark, missing values are "."

#create the encounter histories
trans <- t(dethist) #first transpose so that data for a single site falls in one column instead of one row
ch <- pasteCols(trans, sep="") #now collapse each column into a capture history (ch)
ch <- as.character(ch)

# add covariates - no need to "unlist" with RMark
#on detection

# standardize here! orrrrr not needed, Mark does it for you
jdate <- LIdet[,26:49] 
jdate <- normalize(jdate, method = "range")
jdate[is.na(jdate)] <- 0
tmin = LIdet[,50:73]
tmin <- normalize(tmin, method = "range")
tmin[is.na(tmin)] <- 0
wind <- LIdet[,74:97]
wind <- normalize(wind, method = "range")
wind[is.na(wind)] <- 0
#precip <- LIdet[,98:121]
#precip <- normalize(precip, method = "range")
#precip[is.na(precip)] <- 0
east = LIdet[,122:145]
east <- normalize(east, method = "range")
east[is.na(east)] <- 0

library(corrplot)
pcovs = data.frame(jdate, tmin, wind, east)
C = cor(pcovs, use = "pairwise.complete.obs")
corrplot(C, method = 'number')
corrplot.mixed(C, lower = "number", upper = "circle")


#on to theta
#theta covs -
dist_wet <- sitecovs[,2:4]
dist_wet <-normalize(dist_wet, method = "range")
dist_wet[is.na(dist_wet)] <- 0
dist_fresh <- sitecovs[,5:7]
dist_fresh <-normalize(dist_fresh, method = "range")
dist_fresh[is.na(dist_fresh)] <- 0
dist_edge = sitecovs[,8:10]
dist_edge = normalize(dist_edge, method = "range")
dist_edge[is.na(dist_edge)] <- 0
#dist_rd = sitecovs[,11:13] #needs edits
#dist_rd = normalize(dist_rd, method = "range")
#dist_rd[is.na(dist_rd)] <- 0
dist_shore = sitecovs[,50:52]
dist_shore = normalize(dist_shore, method = "range")
dist_shore[is.na(dist_shore)] <- 0
forest = sitecovs[,14:16]
#forest = normalize(forest, method = "range")
forest[is.na(forest)] <- 0
wetland = sitecovs[,17:19]
#wetland = normalize(wetland, method = "range")
wetland[is.na(wetland)] <- 0
#water = sitecovs[,20:22]           needs edits
#water = normalize(water, method = "range")
#water[is.na(water)] <- 0
#ag = sitecovs[,23:25]              needs edits
#ag = normalize(ag, method = "range")
#ag[is.na(ag)] <- 0
dev = sitecovs[,26:28]
#dev = normalize(dev, method = "range")
dev[is.na(dev)] <- 0
grass = sitecovs[,29:31]
#grass = normalize(grass, method = "range")
grass[is.na(grass)] <- 0
devlow = sitecovs[,32:34]
#devlow = normalize(devlow, method = "range")
devlow[is.na(devlow)] <- 0
devmed = sitecovs[,35:37]
#devmed = normalize(devmed, method = "range")
devmed[is.na(devmed)] <- 0
devhigh = sitecovs[,38:40]
#devhigh = normalize(devhigh, method = "range")
devhigh[is.na(devhigh)] <- 0
devopen = sitecovs[,41:43]
#devopen = normalize(devopen, method = "range")
devopen[is.na(devopen)] <- 0
psm2 = sitecovs[,44:46]
psm2 = normalize(psm2, method = "range")
psm2[is.na(psm2)] <- 0


# on psi
#unit covs - standardize here!
CClvl = LIdet$CC_Lvl
CClvl = as.factor(CClvl)
CC = LIdet$CC
CC <- normalize(CC, method = "range")
MPS = LIdet$MPS
#MPS <- normalize(MPS, method = "range")
NND = LIdet$Patch_NND
#NND <- normalize(NND, method = "range")
EdgeD <- LIdet$EdgeDensity
PA <- LIdet$Paratio
#landscape is now in proportions so 0 to 1
ForestUnit = LIdet$Forest
#ForestUnit <- normalize(ForestUnit, method = "standardize")
AgUnit = LIdet$Ag
#AgUnit <- normalize(AgUnit, method = "standardize")
ImpUnit = LIdet$Imp
#ImpUnit <- normalize(ImpUnit, method = "standardize")
# GrassUnit = LIdet$Grass
# GrassUnit <- normalize(GrassUnit, method = "standardize")
# ShrubUnit = LIdet$Shrub
# ShrubUnit <- normalize(ShrubUnit, method = "standardize")
HerbUnit = LIdet$Grass.Shrub
#HerbUnit <- normalize(HerbUnit, method = "standardize")
NWIfresh = LIdet$NWI_Fresh
#NWIfresh <- normalize(NWIfresh, method = "standardize")
NWIall = LIdet$NWI_All
#NWIall <- normalize(NWIall, method = "standardize")



pscovs = data.frame(CC,MPS,NND,EdgeD,PA,ForestUnit,AgUnit,ImpUnit,HerbUnit,NWIfresh,NWIall)
str(pscovs)
C = cor(pscovs)
corrplot(C, method = 'number')
corrplot.mixed(C, lower = "number", upper = "circle")
#MPS and CC, Edge D and CC, NWI all with fresh


#this next part formats for Program Mark/RMark where each line must end with a ";"
#and the file must be .inp
#can add in covariates if you have them
#create capture history
capture.hist <- NULL
for (q in 1:length(ch)){
  capture.hist[q] <- cbind(paste(c(
    "/*", LIdet[q,1],"*/", # this is the site name and the /* */ are comments
    ch[q],"  ", 1, " ",     # this is the capture history and the frequency of that history (1)
    #CClvl[q],"  ",            #site covs
    CC[q],"  ", 
    MPS[q],"  ", #replace these with other Psi covariate names 
    NND[q],"  ",
    PA[q],"  ",
    EdgeD[q],"  ",
    ForestUnit[q],"  ", 
    AgUnit[q],"  ", 
    ImpUnit[q],"  ", 
    HerbUnit[q],"  ",
    NWIfresh[q],"  ",
    NWIall[q],"  ", 
    paste(jdate[q,],collapse = " "),"  ", #detection covs
    paste(tmin[q,],collapse = " "),"  ",
    paste(wind[q,],collapse = " "),"  ",
    #paste(precip[q,],collapse = " "),"  ",
    paste(east[q,],collapse = " "),"  ",
    paste(dist_wet[q,],collapse = " "),"  ",
    paste(dist_fresh[q,],collapse = " "),"  ",#theta covs
    paste(dist_edge[q,],collapse = " "),"  ",
    #paste(dist_rd[q,],collapse = " "),"  ",
    paste(dist_shore[q,],collapse = " "),"  ",
    paste(psm2[q,],collapse = " "),"  ", 
    paste(forest[q,],collapse = " "),"  ",
    paste(wetland[q,],collapse = " "),"  ",
    #paste(ag[q,],collapse = " "),"  ",
    paste(dev[q,],collapse = " "),"  ",
    paste(devlow[q,],collapse = " "),"  ",
    paste(devmed[q,],collapse = " "),"  ",
    paste(devhigh[q,],collapse = " "),"  ",
    paste(devopen[q,],collapse = " "),"  ",
    paste(grass[q,],collapse = " "),"  ",
    ";"),collapse=""))
}
capture.hist

#write a .inp output file that I can read into MARK/RMark
write(capture.hist, file="simple_bat_hist.inp")


