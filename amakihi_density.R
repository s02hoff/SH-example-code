#
# Estimating Density of Amakihi -------------------------------------------
#


# Libraries & Setup -------------------------------------------------------

rm(list = ls())
library(Distance)


# Data Import and Exploration ---------------------------------------------

honeycreep <- 
  read.csv('data/amakihi.csv')

head(honeycreep)
summary(honeycreep)
str(honeycreep)

names(honeycreep) <-
  c("Study Area", "Region.Label",
    "Area","Sample.Label", "Effort",
    "distance", "OBS", "MAS", "HAS")

#Make HAS a factor
honeycreep$HAS <- 
  as.factor(honeycreep$HAS)

#rescale MAS so it is comparable with the other covariates
honeycreep$MAS <- honeycreep$MAS / sd(honeycreep$MAS, na.rm = TRUE)


max(honeycreep$distance, na.rm = T)
#250m

par(mfrow=c(1,3))

hist(honeycreep$distance,
     main="Default Binning",
     xlab = "Distance")

hist(honeycreep$distance,
     breaks = seq(0,250, length.out = 20),
     main="20 equal size bins",
     xlab = "Distance")

hist(honeycreep$distance,
     breaks = seq(0,250, length.out = 25),
     main="25 equal size bins",
     xlab = "Distance")

hist(honeycreep$distance,
     breaks = c(0,260),
     main="Custom Bins",
     xlab = "Distance")



par(mfrow=c(1,2))

boxplot(
  formula = honeycreep$distance~honeycreep$OBS, 
  xlab="Observer",
  ylab="Distance(m)")

boxplot(
  formula = honeycreep$distance~honeycreep$Region.Label, 
  xlab="Survey Period", 
  ylab="Distance(m)")

boxplot(
  formula = honeycreep$distance~honeycreep$HAS, 
  xlab="HAS", 
  ylab="Distance(m)")

plot(honeycreep$distance~honeycreep$MAS,
     xlab="MAS",
     ylab="Distance(m)")



# Data Analysis - single pooled function ----------------------------------

#run 4 typical models without covariates

hc.hn.herm <- 
  ds(honeycreep, 
     truncation=150, 
     transect="point",
     key="hn", 
     adjustment="herm",
     convert.units = 0.01)

hc.hn.cos <- 
  ds(honeycreep, 
     truncation=150, 
     transect="point",
     key="hn", 
     adjustment="cos", 
     convert.units = 0.01)

hc.uni.cos <- 
  ds(honeycreep, 
     truncation=150, 
     transect="point",
     key="unif", 
     adjustment="cos", 
     convert.units = 0.01)

hc.haz.simp <- 
  ds(honeycreep, 
     truncation=150, 
     transect="point",
     key="hr", 
     adjustment="poly", 
     convert.units = 0.01)


summary1 <- 
  summarize_ds_models(
    hc.haz.simp,
    hc.uni.cos,
    hc.hn.cos,
    hc.hn.herm,
    output = "plain")


View(summary1)
#half-normal with cosine adjustment fits best 

par(mfrow=c(1,3))
plot(hc.hn.cos,
     main="amakihi data pooled, \nhalf-normal cosine detection function", 
     pdf = T)


plot(hc.haz.simp,
     main="amakihi data pooled, \nhazard-rate simple polynomial detection function",
     breaks = c(0.0, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 42.5), 
     pdf = T)

plot(hc.uni.cos,
     main="amikihi data pooled, \nuniform cosine detection function", 
     breaks = c(0.0, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 42.5),
     pdf = T)



par(mfrow=c(1,3))

fit.hc.hn.cos <- 
  ddf.gof(hc.hn.cos$ddf)

fit.hc.uni.cos <- 
  ddf.gof(hc.uni.cos$ddf)

fit.hc.haz.simp <- 
  ddf.gof(hc.haz.simp$ddf)

fit.hc.hn.cos$dsgof$CvM

fit.hc.uni.cos$dsgof$CvM

fit.hc.haz.simp$dsgof$CvM

hc.hn.cos$dht$individuals$summary
hc.hn.cos$dht$individuals$D



# Observer as covariate ---------------------------------------------------

hc.haz.obs <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hr",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ OBS
  )

hc.hn.obs <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hn",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ OBS
  )


summary2 <- 
  summarize_ds_models(
    hc.hn.obs,
    hc.haz.obs,
    hc.hn.cos,
    output = "plain")
View(summary2)

#Observer doesnt seem to be important, basic HN model is still on top



# MAS as covariate --------------------------------------------------------


hc.haz.mas <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hr",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ MAS
  )

hc.hn.mas <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hn",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ MAS
  )

summary3 <- 
  summarize_ds_models(
    hc.hn.mas,
    hc.haz.mas,
    hc.hn.cos,
    output = "plain")
View(summary3)

#HN with MAS is now on top but barely better than base model



# HAS as covariate --------------------------------------------------------


hc.haz.has <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hr",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ HAS
  )

hc.hn.has <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hn",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ HAS
  )

summary4 <- 
  summarize_ds_models(
    hc.hn.has,
    hc.haz.has,
    hc.hn.cos,
    hc.hn.mas,
    output = "plain")
View(summary4)
#HN with MAS model still best 



# Interaction covariates --------------------------------------------------

hc.haz.has.obs <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hr",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ HAS + OBS
  )

hc.hn.has.obs <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hn",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ HAS + OBS
  )

hc.haz.mas.obs <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hr",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ MAS + OBS
  )

hc.hn.mas.obs <-
  ds(
    honeycreep,
    truncation = 150,
    transect = "point",
    key = "hn",
    adjustment = "poly",
    convert.units = 0.01,
    formula = ~ MAS + OBS
  )


summary5 <- 
  summarize_ds_models(
    hc.hn.has.obs,
    hc.haz.has.obs,
    hc.haz.mas.obs,
    hc.hn.mas.obs,
    hc.hn.cos,
    hc.hn.mas,
    output = "plain")
View(summary5)

#HN MAS + OBS is now top model, not far behind is HN HAS + OBS

#summarize all models together

summary6 <- 
  summarize_ds_models(
    hc.hn.has.obs,
    hc.haz.has.obs,
    hc.haz.mas.obs,
    hc.hn.mas.obs,
    hc.hn.cos,
    hc.hn.mas,
    hc.haz.mas,
    hc.hn.obs,
    hc.haz.obs,
    hc.haz.has,
    hc.hn.has,
    output = "plain")
View(summary6)
#HN MAS + OBS stays top model when comparing all models
#Density estimate is now lower at 5.7/ha (CI = 5.32-6.11, CV = 4%)

#Plots for different OBS
par(mfrow=c(1,3))

plot(hc.hn.mas.obs,
     main = "Amikihi Data, \ndetection function for Observer TJS",
     subset = OBS == "TJS",
     ylim = c(0, 1.2))

plot(hc.hn.mas.obs,
     main = "Amikihi, \ndetection function for Observer SGF",
     subset = OBS == "SGF",
     ylim = c(0, 1.2))

plot(hc.hn.mas.obs,
     main = "Amikihi, \ndetection function for Observer TKP",
     subset = OBS == "TKP",
     ylim = c(0, 1.2))

plot(hc.hn.mas.obs,
     main = "Amikihi, \ndetection function",
     ylim = c(0, 1.2))



par(mfrow=c(1,1))
covar.fit <-
  ddf.gof(hc.hn.mas.obs$ddf)
message <-
  paste(
    "Cramer von-Mises W=",
    round(covar.fit$dsgof$CvM$W, 3),
    "\nP=",
    round(covar.fit$dsgof$CvM$p, 3)
  )

text(0.6, 0.1, message, cex=0.8)

hc.hn.mas.obs$dht$individuals$summary
hc.hn.mas.obs$dht$individuals$D






