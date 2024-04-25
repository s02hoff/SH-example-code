#######################
#Running some RMark models#
#######################
rm(list=ls())

#load libraries
library(RMark) # you also need to have MARK installed on your machine


#now to RMark
bat.convert <- convert.inp("simple_bat_hist.inp", 
                            covariates = c("CC", "MPS","EdgeD", "ForestUnit", "AgUnit", "ImpUnit", "HerbUnit", "NWIfresh","date11", "date12", "date13", "date14", "date15", "date16", "date17","date18", "date21", "date22", "date23", "date24", "date25", "date26","date27", "date28", "date31", "date32", "date33", "date34", "date35","date36", "date37", "date38", "tmin11", "tmin12", "tmin13", "tmin14", "tmin15", "tmin16", "tmin17","tmin18", "tmin21", "tmin22", "tmin23", "tmin24", "tmin25", "tmin26","tmin27", "tmin28", "tmin31", "tmin32", "tmin33", "tmin34", "tmin35","tmin36", "tmin37", "tmin38", "wind11", "wind12", "wind13", "wind14", "wind15", "wind16", "wind17","wind18", "wind21", "wind22", "wind23", "wind24", "wind25", "wind26","wind27", "wind28", "wind31", "wind32", "wind33", "wind34", "wind35","wind36", "wind37", "wind38","dettype11","dettype12","dettype13","dettype14","dettype15","dettype16","dettype17","dettype18","dettype21","dettype22","dettype23","dettype24","dettype25","dettype26","dettype27","dettype28","dettype31","dettype32","dettype33","dettype34","dettype35","dettype36","dettype37","dettype38","date211", "date212", "date213", "date214", "date215", "date216", "date217","date218", "date221", "date222", "date223", "date224", "date225", "date226","date227", "date228", "date231", "date232", "date233", "date234", "date235","date236", "date237", "date238","dwet1", "dwet2", "dwet3", "dedge1", "dedge2", "dedge3", "psm1", "psm2", "psm3", "forest1", "forest2", "forest3", "wet1", "wet2", "wet3","dev1", "dev2", "dev3", "devl1", "devl2", "devl3", "devm1", "devm2", "devm3", "devh1", "devh2", "devh3", "devo1", "devo2", "devo3"))

bat.proc <- process.data(bat.convert, model="MultScalOcc",  mixtures = 8) # mixtures is the number of secondary samples, #groups =c("CClvl")
bat.ddl <- make.design.data(bat.proc) 
str(bat.ddl)


# 
# 
# Calculate probability of detection per sampling night to get effort needed to determine likelihood of presence
p_method <- c(.51,.51, .51)       #  detection prob for each method
p_survey <- 1 - prod(1-p_method)  #  prob detection in each survey
p_overall <- 1 - (1-p_survey)^8   #  prob detection in at least 1 of 8 surveys
cat('p_overall:',p_overall,'\n')
# 
# 
maxsurveys=4; poverall=rep(0, maxsurveys)
pmeth = c(.52, .51, .48)
for (nsurveys in 1:maxsurveys) {
   psurvey = 1 - prod(1-pmeth)
   poverall[nsurveys] = 1 - (1 - psurvey)^nsurveys
 }
 
df=data.frame(surveys=1:maxsurveys, detectionRate=poverall)
ggplot(df, aes(x=surveys, y=detectionRate)) + geom_point() +
   geom_line() + geom_hline(yintercept=0.95) + ylim(0,1)






#look at detection rate varying by time, holding other parameters fixed 

m1 <- mark(bat.proc,bat.ddl, model.parameters = list(
  Psi=list(formula=~1), Theta=list(formula=~1), p=list(formula=~Time)
))
summary(m1)
m1$results$real
#avg = 0.56




timeval <- seq(1,8, length.out=8)
PIMS(m1,"p",simplified=FALSE) 
as.numeric(row.names(bat.ddl$p))
p.predict <- covariate.predictions(m1, indices = c(bat.ddl$p$model.index))
p.predict$estimates
mean(p.predict$estimates$estimate)#0.56
LItime <- data.frame(p.predict$estimates)
write.csv(LItime, "C:\\Users\\sxhoff\\OneDrive - New York State Office of Information Technology Services\\Documents\\R\\Occu\\LC\\LItime.csv", row.names=TRUE)





#Test det covs first
#create function to build models

do.Species<-function()
{
  Psi.dot = list(formula = ~1)
  #Psi.Imp = list(formula = ~ImpUnit)
  Psi.MPS = list(formula = ~MPS)
  #Psi.For = list(formula = ~ForestUnit)
  #Psi.NWI = list(formula = ~NWIfresh)
  Psi.CC = list(formula = ~ CC)
  #Psi.Edge = list(formula = ~ EdgeD)
  #Psi.Ag = list(formula = ~ AgUnit)
  #Psi.Herb = list(formula = ~ HerbUnit)
  
  p.dot = list(formula = ~1)
  p.date = list(formula = ~date)
  #p.tmin = list(formula = ~tmin)
  p.wind = list(formula = ~wind)
  #p.type = list(formula = ~dettype)
  #p.date2 = list(formula = ~date2)
  #p.time = list(formula = ~Time)
  #p.1 = list(formula = ~date + tmin)
  #p.2 = list(formula = ~date + wind)
  #p.3 = list(formula = ~date + dettype)
  #p.4 = list(formula = ~tmin + wind)
  #p.5 = list(formula = ~tmin + dettype)
  #p.6 = list(formula = ~wind + dettype)
  #p.7 = list(formula = ~tmin + date2)
  #p.8 = list(formula = ~wind + date2)
  #p.9 = list(formula = ~dettype + date2)
  
  #Theta.dot = list(formula = ~1)
  #Theta.psm = list(formula = ~psm)
  #Theta.f = list(formula = ~forest)
  #Theta.dwet = list(formula = ~dwet)
  #Theta.dedge = list(formula = ~dedge)
  Theta.dev = list(formula = ~dev)
  #Theta.devl = list(formula = ~devl)
  #Theta.devm = list(formula = ~devm)
  #Theta.devh = list(formula = ~devh)
  #Theta.devo = list(formula = ~devo)

  
  occ=create.model.list("MultScalOcc")
  return(mark.wrapper(occ,data=bat.proc,ddl=bat.ddl, adjust=FALSE,realvcv=TRUE, use.initial = TRUE))
}
#options = "SIMANNEAL"
bat.results <- do.Species()

# Output model table and estimates
bat.results$model.table
bat.results[[as.numeric(rownames(bat.results$model.table[4,]))]]$results$real
bat.results[[as.numeric(rownames(bat.results$model.table[6,]))]]$results$beta


#tested first with p covs
#holding Psi with Forest and theta with forest
#tmin and date in top mods 


#theta covs
#dev is on top

#tested psi covs last
#no relationships

top <- bat.results$Psi.dot.Theta.f.p.dot
summary(top)
summary.mark(top, se=TRUE)
coef(top)

fc <- find.covariates(top,bat.convert)
design <- fill.covariates(top,fc)
cr <- compute.real(top, design=design)

date.values <- seq(0,1, length.out=100)
covariate.predictions(bat.results, data=data.frame(date.values), indices=c(bat.ddl$p$model.index))

p.estimates = model.average(bat.results,"p", vcv = TRUE, drop = TRUE) # drop = true will drop any models with non-positive variance for betas
print(unique(p.estimates$estimate))
mean(p.estimates$estimate$estimate) #0.6655219


psi.estimates = model.average(bat.results,"Psi",vcv = TRUE, drop = TRUE)
print(unique(psi.estimates$estimates$estimate)) 

theta.estimates = model.average(bat.results,"Theta", vcv = TRUE, drop = TRUE)
print(unique(theta.estimates$estimate))
mean(theta.estimates$estimate$estimate)#0.6005004

#Predict p by survey night from model averages
p.predict <- covariate.predictions(bat.results, indices = c(bat.ddl$p$model.index))
p.predict$estimates
LItime1 <- data.frame(p.predict$estimates)
write.csv(LItime1, "C:\\Users\\sxhoff\\OneDrive - New York State Office of Information Technology Services\\Documents\\R\\Occu\\LC\\LItime.csv", row.names=TRUE)


#PLOTS FOR P-----

#Easting
#min and max 631640,762225,
date.values <- seq(0,1, length.out=100)
PIMS(bat.results[[1]],"p",simplified=FALSE) 
as.numeric(row.names(bat.ddl$p))
e.predict <- covariate.predictions(bat.results, data=data.frame(date12=date.values), indices=c(6))#indices=c(bat.ddl$p$model.index)
e.predict$estimates
259900

plot(e.predict$estimates$covdata,e.predict$estimates$estimate,type='l',ylim=c(0,1))
lines(east.values,e.predict$estimates$lcl,pch=16,col='red')
lines(east.values,e.predict$estimates$ucl,pch=16,col='red')
d2 <- seq(631350,762230, length.out=100)
occu_east_pred_df <- data.frame(Predicted = e.predict$estimates$estimate,
                                lower = e.predict$estimates$lcl,
                                upper = e.predict$estimates$ucl,
                                data=data.frame(east11=d2),isl="Long Island")

# Plot with ggplot
occu_east_pred_plot <- ggplot(occu_east_pred_df, aes(x = east11, y = Predicted, fill=isl)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "solid") +
  geom_path(size = 1) +
  labs(x = "Easting (UTM)", y = "Detection probability") + theme(legend.title = element_blank(),legend.box.background = element_rect(color = "black"), legend.background = element_blank(),legend.position = c(0.60, 0.85), legend.justification = c(1,0)) +
  theme_classic() + coord_cartesian(ylim = c(0,1)) +
  theme(legend.position = "top",text = element_text(size = 25))  + scale_fill_manual(values=c("black"))
occu_east_pred_plot


#+ ggtitle("LI") + theme(plot.title=element_text(margin=margin(t=40,b=-30)))
occu_east_pred_plot


#PLOTS for THETA-----
#Dev on Theta
dev.values <- seq(0,1,length.out=66)
PIMS(bat.results[[1]],"Theta",simplified=FALSE) 
as.numeric(row.names(bat.ddl$Theta))
dev.predict <- covariate.predictions(bat.results, data=data.frame(dev2=dev.values), indices = c(3))#indices=c(bat.ddl$p$model.index)
dev.predict$estimates
#devvals <- seq(0,30, length.out=100)#for acres
occu_dev_pred_df <- data.frame(Predicted = dev.predict$estimates$estimate,
                               lower = dev.predict$estimates$lcl,
                               upper = dev.predict$estimates$ucl,
                               data=data.frame(dev2=dev.values),occ=liocc)

#dont round
#occu_dev_pred_df <- occu_dev_pred_df %>% 
#  mutate_if(is.numeric, round, digits=2)

#add raw data
lires <- read.csv("LISiteResult.csv", header = TRUE, na.strings = 'NA')
str(lires)
lidev <- lires$Dev
liocc <- lires$Occu
ggplot() + geom_point(lires, mapping = aes(x=Dev, y = Occu))

#plot
occu_dev_pred_plot <- ggplot(occu_dev_pred_df, aes(x = dev2, y = Predicted))  + geom_point(mapping = aes(x=dev2, y=occ)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "solid") + geom_path(size = 1) +
  labs(x = "Proportion of development within 200 m", y = "Local occupancy probability") +
  theme_classic() + coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(size = 20))
occu_dev_pred_plot
#bellisimo!


#forest on theta
f.values <- seq(0,1,length.out=100)
f.predict <- covariate.predictions(bat.results, data=data.frame(forest1=f.values), indices = c(2))#indices=c(bat.ddl$p$model.index)
f.predict$estimates

occu_for_pred_df <- data.frame(Predicted = f.predict$estimates$estimate,
                               lower = f.predict$estimates$lcl,
                               upper = f.predict$estimates$ucl,
                               data=data.frame(forest=f.values))

occu_for_pred_plot <- ggplot(occu_for_pred_df, aes(x = forest, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "solid") +
  geom_path(size = 1) +
  labs(x = "Proportion of forest within 200 m", y = "Local occupancy probability") +
  theme_classic() + coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(size = 20))
occu_for_pred_plot


#Dev open on theta
devo.values <- seq(0,1,length.out=100)
devo.predict <- covariate.predictions(bat.results, data=data.frame(devl1=devo.values), indices = c(2))#indices=c(bat.ddl$p$model.index)
devo.predict$estimates

occu_devo_pred_df <- data.frame(Predicted = devo.predict$estimates$estimate,
                               lower = devo.predict$estimates$lcl,
                               upper = devo.predict$estimates$ucl,
                               data=data.frame(devo=devo.values))

occu_devo_pred_plot <- ggplot(occu_devo_pred_df, aes(x = devo, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "solid") +
  geom_path(size = 1) +
  labs(x = "Proportion of open development within 200 m", y = "Local occupancy probability") +
  theme_classic() + coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(size = 20))
occu_devo_pred_plot

#Dev low on Theta
devl.values <- seq(0,1,length.out=100)
devl.predict <- covariate.predictions(bat.results, data=data.frame(devl3=devl.values), indices = c(4))#indices=c(bat.ddl$p$model.index)
devl.predict$estimates

plot(devl.predict$estimates$covdata,devl.predict$estimates$estimate,type='l',ylim=c(0,1))
lines(devl.values,devl.predict$estimates$lcl,pch=16,col='red')
lines(devl.values,devl.predict$estimates$ucl,pch=16,col='red')

#Dev high on theta
devh.values <- seq(0,1,length.out=100)
devh.predict <- covariate.predictions(bat.results, data=data.frame(devh2=devh.values), indices = c(3))#indices=c(bat.ddl$p$model.index)
devh.predict$estimates




#PLOTS for PSI-----
#NONE =(

#MPS for Psi - nope
imp.values <- seq(0,1, length.out=100)
PIMS(bat.results[[1]],"Psi",simplified=FALSE) 
as.numeric(row.names(bat.ddl$Psi))
imp.predict <- covariate.predictions(bat.results, data=data.frame(ImpUnit=imp.values), indices = c(1))#indices=c(bat.ddl$p$model.index)
imp.predict$estimates



  # this creates a lot of temporary files in your working directory, so before I'm done for the day I usually run the following to remove them
rm(list=ls())
cleanup(ask=FALSE)


