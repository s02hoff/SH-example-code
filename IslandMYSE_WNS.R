#ISLAND MYSE WNS-----


#load packages
library(lubridate)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(multcomp)
library(tidyr)
library(lme4)
library(dplyr)
library(ggpubr)
library(scales)
library(effects)
library(sjPlot)
library(car)
library(rcompanion)
library(scales)
library(AICcmodavg)
library(knitr)
library(glmmTMB)
library(DHARMa)

#load data file
setwd("/Users/samso/Dropbox/R Files")

wnsdat <- read.csv("PdData_Island_Mainland_Combo.csv", header = TRUE, na.strings = 'NA')
options(width=100)
print(head(wnsdat))
str(wnsdat)

#summarise data 
group_by(wnsdat, Location) %>%
  dplyr::summarize(n=n()) #145 LI records, 530 mainland
group_by(wnsdat, Location, Season) %>%
  dplyr::summarize(n=n())
#By Season
#Island: 68 fall, 57 spring, 20 summer
#Main: 355 fall, 175 spring
with(wnsdat, table(Location,SurveySite))

#split by disease stage (invasion=0-2, post-invasion=3+)
group_by(wnsdat, State,Site,Season,Group) %>%
  dplyr::summarize(n=n())


#WNS Models-----
##Prevalence-----

#subset to fall and spring data, remove summer 
wnsdat2 <- subset(wnsdat, Season %in% c('Fall','Spring'))%>% droplevels()#655 obs
str(wnsdat2)

#summarise prevalence data by site, season and group for appendix table
Prevsum<- wnsdat2 %>% 
  group_by(Prev,Site,Season,Group) %>%
  dplyr::summarise(n = n()) 


###Island prevalence-----
#Figure 2a
#subset dataframe to get just the island data 
isldat <- subset(wnsdat2, Location %in% c('Island'))%>% droplevels()#125 obs
str(isldat)

#summarise data
with(isldat, table(Season, Prev))
#Positives = 21% overall, fall = 3/68 (0.04%); spring = 23/57 (40%)
isldat %>% 
  group_by(Prev, Status, Season) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

#model with time (pdate) and site (island) as fixed effects, site as a random effect
P.Isl = glmmTMB(Prev ~ Pdate + Site + (1|SurveySite), data=isldat, family="binomial")
summary(P.Isl)
plot_model(P.Isl, type = "pred", terms = c("Pdate[all]","Site"))
tab_model(P.Isl, p.style = "numeric_stars", transform = NULL)


#model predictions
newdata1 <- crossing(Pdate =seq(10,18,length.out=100), Site = c("LI","MV","N"), SurveySite=NA)
newdata1
newdata1$prev <- predict(P.Isl, newdata = newdata1, type = "response")
newdata1
newdata2 <- cbind(newdata1, predict(P.Isl, newdata = newdata1, type = "link", se = TRUE))
newdata2 <- within(newdata2, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
head(newdata2)

#newdata2 <- newdata2 %>% 
#  mutate_if(is.numeric, round, digits=2)

o3 <- isldat %>%
  group_by(Pdate,Site) %>%
  summarise(prev=mean(Prev,na.rm=T),sd=sd(Prev,na.rm=T),n=n()) %>%
  mutate(se = (prev*(1-prev)/n)^.5)


#Plot Prev by Date & Site
#if you want proportion of positive by day, use PdateProp
ggplot() + geom_point(data=o3, aes(Pdate,prev, color=Site, size=n),alpha=.5) + geom_ribbon(data = newdata2,aes(x=Pdate, y=PredictedProb,ymin = LL,ymax = UL, fill = Site), alpha = 0.2) + geom_line(data=newdata2,aes(x=Pdate,y=PredictedProb,colour = Site),size = 1) + scale_fill_manual(values=c( "#009E73", "#F0E442", "#D55E00")) + theme_classic() + scale_color_manual(values=c( "#009E73", "#F0E442", "#D55E00")) + theme(text = element_text(size = 20), legend.text=element_text(size=15),legend.position = "top") + labs(y=expression(~italic("P. destructans")~"prevalence"))+ scale_x_continuous(breaks = c(10.00000,12.02020,14.04040,16.06061,18.00000),labels = c("10","12", "2", "4", "6"))+ scale_size(guide = "none") + theme(legend.title=element_blank()) + xlab("month")

#axis.title.x=element_blank()
###Disease phase-----
#Figure 2b

P.Stage = glmmTMB(Prev ~ Pdate+Group*Location + (1|SurveySite), data=wnsdat2, family="binomial")
summary(P.Stage)
plot_model(P.Stage, type = "pred", terms = c("Pdate[all]","Group","Location"))
pstage_simresid <- simulateResiduals(P.Stage)
plot(pstage_simresid)

model_means_prev = emmeans(P.Stage, ~ Group*Location)
as.data.frame(model_means_prev)

model_means_cld <- cld(object = model_means_prev,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)
#post-invasion island not sig different from mainland invasion years
pairs(model_means_prev)


#get model predictions
nd <- crossing(Pdate =seq(8,18,length.out=100), Location = c("Island","Main"), Group = c("Invasion","Post-invasion"), SurveySite = NA)

nd2 <- cbind(nd, predict(P.Stage, newdata = nd, type = "link", se = TRUE))
nd2 <- within(nd2, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
head(nd2)

#get raw data points for figure
wnsprev <- read.csv("PdData_Island_Mainland_Combo.csv", header = TRUE, na.strings = 'NA') #load og data, keep prev numeric
wnsprev <- subset(wnsprev, Season %in% c('Fall','Spring'))%>% droplevels()
pd1 <- wnsprev %>%
  group_by(Pdate,Location,Group, YSW) %>%
  summarise(prev=mean(Prev,na.rm=T),sd=sd(Prev,na.rm=T),n=n()) %>%
  #  summarise_at(vars(prev),funs(prev=mean(na.rm=T),sd(na.rm=T),n())) %>%
  mutate(se = (prev*(1-prev)/n)^.5)

pd1<- pd1 %>% drop_na(Group)
str(pd1)


nd2 %>%
  mutate(across(Group, factor, levels=c("Invasion","Post-invasion")))

#plot one - invasion years
pd1inv <- subset(pd1, Group %in% c('Invasion'))%>% droplevels()
pd1inv$YSW <- as.factor(pd1inv$YSW)
nd2inv <- subset(nd2, Group %in% c('Invasion'))%>% droplevels()

invplot <- ggplot () + geom_ribbon(data=nd2inv, aes(x=Pdate, y=PredictedProb,ymin = LL,ymax = UL, fill = Location), alpha = 0.2) + geom_line(data=nd2inv,aes(x=Pdate,y=PredictedProb,color = Location),size = 1) + geom_point(data=pd1inv, aes(x=Pdate,y=prev, size = n, alpha=YSW, color = Location), position = "jitter") + scale_size(guide = "none") + theme_classic() + scale_color_manual(values=c("deepskyblue4","deeppink4"))+ scale_fill_manual(values=c("deepskyblue4","deeppink4")) + labs(y=expression(~italic("P. destructans")~"prevalence"),x="month") + scale_x_continuous(breaks = c(8.000000,10.000000,12.02020,14.04040,16.06061,18.00000),labels = c("8","10","12", "2","4","6")) + theme (legend.position = "none", text = element_text(size = 20)) + guides(fill = "none", color="none") + scale_alpha_discrete(range = c(0.1, 0.3))

invplot

#plot two - post-invasion years
pd1post <- subset(pd1, Group %in% c('Post-invasion'))%>% droplevels()
pd1post$YSW <- as.factor(pd1post$YSW)
nd2post <- subset(nd2, Group %in% c('Post-invasion'))%>% droplevels()

postplot <- ggplot () + geom_ribbon(data=nd2post, aes(x=Pdate, y=PredictedProb,ymin = LL,ymax = UL, fill = Location), alpha = 0.2) + geom_line(data=nd2post,aes(x=Pdate,y=PredictedProb,color = Location),size = 1) + geom_point(data=pd1post, aes(x=Pdate,y=prev, size = n, alpha=YSW, color = Location), position = "jitter") + scale_size(guide = "none") + theme_classic() + scale_color_manual(values=c("deepskyblue4","deeppink4"))+ scale_fill_manual(values=c("deepskyblue4","deeppink4")) + scale_x_continuous(breaks = c(8.000000,10.000000,12.02020,14.04040,16.06061,18.00000),labels = c("8","10","12", "2","4","6")) + theme (legend.position = "none", text = element_text(size = 20), axis.title.y = element_blank(), axis.text.y = element_blank()) + xlab("month")+ guides(fill = "none", color="none")+ scale_alpha_discrete(range = c(0.4, 1.0))

postplot

invplot + postplot

#plot with YSW legend
pd1$YSW <- as.factor(pd1$YSW)
p.all.ysw <- ggplot () + facet_wrap(~Group) + geom_ribbon(data=nd2, aes(x=Pdate, y=PredictedProb,ymin = LL,ymax = UL, fill = Location), alpha = 0.2) + geom_line(data=nd2,aes(x=Pdate,y=PredictedProb,color = Location),size = 1) + geom_point(data=pd1, aes(x=Pdate,y=prev, size = n, alpha=YSW, color = Location), position = "jitter") + theme_classic() + theme(legend.position = "top",text = element_text(size = 20), legend.text=element_text(size=15), legend.title = element_text(size=15)) + scale_color_manual(values=c("deepskyblue4","deeppink4"))+ scale_fill_manual(values=c("deepskyblue4","deeppink4")) + labs(y=expression(~italic("P. destructans")~"prevalence"),x="month") + scale_x_continuous(breaks = c(8.000000,10.000000,12.02020,14.04040,16.06061,18.00000),labels = c("8","10","12", "2","4","6")) + guides(alpha=guide_legend(nrow=1), fill="none", color="none") + scale_size(guide = "none")

p.all.ysw
 
#plot with location legend
p.all.loc <- ggplot () + facet_wrap(~Group) + geom_ribbon(data=nd2, aes(x=Pdate, y=PredictedProb,ymin = LL,ymax = UL, fill = Location), alpha = 0.2) + geom_line(data=nd2,aes(x=Pdate,y=PredictedProb,color = Location),size = 1) + geom_point(data=pd1, aes(x=Pdate,y=prev, size = n, alpha=YSW, color = Location), position = "jitter") + theme_classic() + theme(legend.position = "top",text = element_text(size = 20), legend.text=element_text(size=15),legend.title=element_blank()) + scale_color_manual(values=c("deepskyblue4","deeppink4"))+ scale_fill_manual(values=c("deepskyblue4","deeppink4")) + labs(y=expression(~italic("P. destructans")~"prevalence"),x="month") + scale_x_continuous(breaks = c(8.000000,10.000000,12.02020,14.04040,16.06061,18.00000),labels = c("8","10","12", "2","4","6")) + guides(alpha="none") + scale_size(guide = "none")

p.all.loc


library(gridExtra)
library(cowplot)
library(grid)
#extract legend from plot
grobs <- ggplotGrob(p.all.loc)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
s1<-plot_grid(p.all.ysw, nrow = 1, labels=c('', ''),rel_heights = c(1/2,1/2))#label_x =.07
s1
s2<-plot_grid(legend,s1, nrow = 2, labels=c('', ''),rel_heights = c(1/9, 8/9))#label_x =.07
s2



####YSW-----
#Figure 2a
#just use spring data when prevalence is at its highest
wnsdat3 <- subset(wnsdat, Season %in% c('Spring'))%>% droplevels()

YSWsum <- wnsdat3 %>% 
  group_by(Prev, YSW, Location) %>%
  dplyr::summarise(n = n()) %>%
  mutate(freq = n / sum(n))

#islands vs mainland over time since WNS arrival

#mixed model with random site effect
P.YSW = glmmTMB(Prev ~ YSW+Location+(1|SurveySite), data=wnsdat3, family="binomial")
summary(P.YSW)
plot(allEffects(P.YSW), type = "response")
plot_model(P.YSW, type = "pred", terms = c("YSW","Location"))
emmeans(P.YSW, pairwise~Location, type="response")
ysw_simresid <- simulateResiduals(P.YSW)
plot(ysw_simresid)

#get model predictions 
newdata1 <- crossing(YSW =seq(0,8,length.out=100), Location = c("Main","Island"),SurveySite= NA)
newdata1
newdata1$prev <- predict(P.YSW, newdata = newdata1, type = "response")
newdata1
newdata2 <- cbind(newdata1, predict(P.YSW, newdata = newdata1, type = "link", se = TRUE))
newdata2 <- within(newdata2, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
head(newdata2)

#Split data to show at what year we are extrapolating for (mainland - no data post YSW 5)
ndmain <- subset(newdata2, Location %in% c('Main'))%>% droplevels()
ndmain1<- subset(ndmain,ndmain$YSW>=	
                   0.00000000 & ndmain$YSW<6.000000)
ndmain2<- subset(ndmain,ndmain$YSW>=6.000000)
ndisl <- subset(newdata2, Location %in% c('Island'))%>% droplevels()

#get data points 
wnsprev <- read.csv("PdData_Island_Mainland_Combo.csv", header = TRUE, na.strings = 'NA') #load original data, keep prev numeric
wnsprev <- subset(wnsprev, Season %in% c('Spring'))%>% droplevels()
pd1 <- wnsprev %>%
  group_by(YSW,Location) %>%
  summarise(prev=mean(Prev,na.rm=T),sd=sd(Prev,na.rm=T),n=n()) %>%
  #  summarise_at(vars(prev),funs(prev=mean(na.rm=T),sd(na.rm=T),n())) %>%
  mutate(se = (prev*(1-prev)/n)^.5)

ggplot() + geom_line(data=ndisl,aes(x=YSW,y=PredictedProb, color = Location)) + geom_ribbon(data=ndisl,aes(x=YSW,ymin = LL,ymax = UL, fill = Location), alpha = 0.2) + geom_line(data=ndmain1,aes(x=YSW, y=PredictedProb,colour = Location),size = 1) + geom_ribbon(data=ndmain1,aes(x=YSW,ymin = LL,ymax = UL, fill = Location), alpha = 0.2)+ geom_line(data=ndmain2,aes(x=YSW, y=PredictedProb,colour = Location),size = 1, linetype = "dashed") + geom_ribbon(data=ndmain2,aes(x=YSW,ymin = LL,ymax = UL, fill = Location), alpha = 0.2) + theme_classic() + theme(text = element_text(size = 20), legend.text=element_text(size=15),legend.position = "top")+ scale_color_manual(values=c("deepskyblue4","deeppink4")) + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + labs(y=expression(~italic("P. destructans")~"prevalence")) + scale_x_continuous(breaks = c(0.00000000,1.00000000,2.00000000,3.00000000,4.00000000,5.00000000,6.00000000,7.00000000,8.00000000),labels = c("0","1","2", "3", "4", "5","6","7","8"))+xlab("Years since WNS arrival") + theme(legend.title=element_blank())+ geom_point(data=pd1, aes(x=YSW,y=prev,color = Location, size = n), alpha=0.5)+ scale_size(guide = "none")



## Infection intensity - fungal loads------

#Remove Nantucket data - only one positive sample
wnsdat4 <- subset(wnsdat2, Site %in% c('LI', 'MV','Main'))%>% droplevels() #576 obs

#Load by Pdate
L.Load = glmmTMB(Log10Load ~ Pdate*Location + (1|Wyear), data=wnsdat4, family = gaussian(link="identity"))
summary(L.Load)
plot_model(L.Load, type = "pred", terms = c("Pdate[all]","Location"))
load_simresid <- simulateResiduals(L.Load)
plot(load_simresid)

#Means by location and date
marginal = lsmeans(L.Load,
                   ~ Pdate*Location)
as.data.frame(marginal)
CLD = cld(marginal,
          alpha   = 0.05,
          Letters = letters,
          adjust  = "sidak")
pairs(marginal,
      alpha=0.05,
      Letters=letters,
      adjust="sidak")

#predict and graph
nd <- crossing(Pdate =seq(8,18,length.out=100), Location = c("Island","Main"), Wyear=NA)
str(nd)
nd2 <- cbind(nd, predict(L.Load, newdata = nd, type = "response", se = TRUE))
nd2 <- within(nd2, {
  PredictedProb <- fit
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)
})
head(nd2)
min(nd2$UL)

nd2$LL[nd2$LL< -6.5] = -6.5
nd2$UL[nd2$UL< -6.5] = -6.5
nd2$fit[nd2$fit< -6.5] = -6.5




ggplot()+ geom_point(data=wnsdat4,aes(x=Pdate,y=Log10Load,color=Location), alpha = 0.5,position="jitter") + geom_ribbon(data=nd2,aes(x=Pdate, y=fit,ymin = LL,ymax = UL, fill = Location), alpha = 0.2) + geom_line(data=nd2,aes(x=Pdate,y=fit,color = Location),size = 1)+ theme_classic()+theme(legend.position = "top",text = element_text(size = 20), legend.text=element_text(size=15), axis.title.x=element_blank(),axis.title.y = element_text(vjust = +2)) + scale_color_manual(values=c("deepskyblue4","deeppink4")) + scale_fill_manual(values=c("deepskyblue4","deeppink4"))+ ylab(expression(Log[10] ~ Pd ~ Loads))+ scale_x_continuous(breaks = c(8.00000,10.00000,12.02020,14.04040,16.06061,18.00000),labels = c("Aug","Oct","Dec","Feb", "April","May")) + theme(legend.title=element_blank())

#+ geom_point(data=wnsdatabs,aes(x=Pdate,y=Log10Load,color=Location),position=position_dodge(width=0.3), alpha = 0.5)





#MYSE decline -----
#import data
epi_dec <- read.csv("epi_dec.csv", header = TRUE, na.strings = 'NA')
options(width=100)
print(head(epi_dec))
str(epi_dec)


#74 sites with at least 2 surveys post WNS
epi_dec1 <- subset(epi_dec, postwns_surv %in% c('2','3','4'))%>% droplevels() 

#summarize by state and site, include minimum decline lambda and winter temps
epidemic = epi_dec1 %>%
  group_by(state,site)%>%
  summarise(min.lambda = min(lambda.post),wtemp=wtemp)

#only remove these if we dont use the filter above 
#remove aunt beck simmons cave from TN
epidemic<-epidemic[!(epidemic$site== "AUNT BECK SIMMONS CAVE"),]
#lets remove the sites from Pa with old data 
epidemic<-epidemic[!(epidemic$site== "EISWERT NO. 2"),]
epidemic<-epidemic[!(epidemic$site== "LAUREL CAVERNS"),]
epidemic<-epidemic[!(epidemic$site== "COON"),]


epi_sum <- epi_dec1 %>%
  group_by(state)%>%
  summarise(min.ysw = min(ysw),max = max(ysw),mean = mean(ysw), n=n())
#YSW range 1-4, mean = 2.51
mean(epi_sum$mean)

epi_sum1 <- epi_dec1 %>%
  group_by(state)%>%
  summarise(min.surv = min(surveydiff),max = max(surveydiff),mean = mean(surveydiff), n=n())
#0-9 years for pre WNS survey


#calculate Lambda mean, sd, and se by state, and average winter temps by state
#winter temp = average temp in C for Dec-Feb, 2011-2020
epi <- epidemic %>%
  group_by(state) %>%
  summarise(mean=mean(min.lambda,na.rm=T),sd=sd(min.lambda,na.rm=T),n=n(),meant=mean(wtemp)) %>%
  mutate(se = (mean*(1-mean)/n)^.5)

#make the mean lambda negative to look at decline since WNS arrival
epi$per <- 1 - epi$mean
epi$meanneg <- -1 * epi$per
epi$ymin <- epi$meanneg - epi$se
epi$ymax <- epi$meanneg + epi$se

#cutoff se at -1.0 for visualization purposes
epi$ymin[epi$ymin > "-1.0"] <- -1.0

#reorder states by average latitude of sites
#MI = 46.51; WI = 43.92; VT = 43.76; NY = 42.82; MA = 42.17; CT = 41.75, IL = 41.34, PA = 40.54; VA = 37.76; TN = 36.29; NC = 35.43
library(forcats)
epi$state <- fct_relevel(epi$state, "MI", "WI", "VT", "NY", "MA", "CT", "IL", "PA", "VA", "TN", "NC")
epidemic$state <- fct_relevel(epidemic$state, "MI", "WI", "VT", "NY", "MA", "CT", "IL", "PA", "VA", "TN", "NC")
#add raw data to the graph
#make lambda negative, and convert to total decline since WNS arrival
epidemic$lambdaneg <- 1 - epidemic$min.lambda
epidemic$lambdaneg <- -1 * epidemic$lambdaneg

p5 <- ggplot(epi, aes(x=state, y=meanneg)) + geom_point(data=epidemic, aes(x=state, y=lambdaneg), alpha=0.5, position=position_jitter(width=0.3),size=2) + 
  geom_errorbar(aes(ymin=ymin, ymax=ymax, color=meant), width=0, size=1) + geom_point(aes(size=n, color=meant), stroke=1) + geom_hline(yintercept=0, linetype="dashed", color = "black") + theme_classic() + scale_size(guide = 'none') + ylab("Percent change after WNS arrival") + theme(axis.title.x = element_blank(), text = element_text(size = 25)) + scale_color_gradient(low = "blue", high = "red", n.breaks=6, labels=c(-8,-5,-2,0,2,5))  
p5 + scale_y_continuous(labels = scales::percent, limits = c(-1.0,0),oob = scales::squish) + theme(legend.position = c(0.95, 0.75), legend.justification = c(1,0), legend.background = element_blank(),legend.text=element_text(size=10),legend.title=element_text(size=15)) + labs(color='Mean winter surface\ntemperature\u00b0C') + guides(color=guide_colorbar(title.position = "top",direction="horizontal",ticks.colour = "NA", frame.colour = "black", frame.linewidth = 0.25)) 


#####Capture probability-----


cap = read.csv("CapProb.csv", header = TRUE, na.strings = 'NA')
options(width=100)
print(head(cap))
cap<- cap %>% drop_na()
cap$Location <- as.factor(cap$Location)
cap$Loc2 <- as.factor(cap$Loc2)
#cap$YSW <- as.factor(cap$YSW)
cap$NN <- cap$Nights * cap$Nets
cap$NNlog <- log10(cap$NN)


cap.glm1 <- glm(Cap ~ YSW * Loc2 + offset(NNlog), data=cap, family = "binomial")
summary(cap.glm1)
cap_simresid1 <- simulateResiduals(cap.glm1)
plot(cap_simresid1)
plot_model(cap.glm1, type = "pred", terms = c("YSW[0:10]","Loc2"))

#model predictions
newdata3 <- crossing(YSW =seq(0,10,length.out=100), Loc2 = c("Mainland","Island"),NNlog=mean(cap$NNlog))
newdata3
newdata3$prev <- predict(cap.glm1, newdata = newdata3, type = "response")
newdata3
newdata4 <- cbind(newdata3, predict(cap.glm1, newdata = newdata3, type = "link", se = TRUE))
newdata4 <- within(newdata4, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
head(newdata4)

rawcap1 <- cap %>%
  group_by(YSW,Loc2,MYSE) %>%
  summarise(prev=mean(Cap,na.rm=T),sd=sd(Cap,na.rm=T),n=n()) %>%
  mutate(se = (prev*(1-prev)/n)^.5)

rawcap1 <- subset(rawcap1, YSW %in% c('0','1','2','3','4','5','6','7','8','9','10'))%>% droplevels()

#Split data to show at what year we are extrapolating for (island - no data post YSW 9)
islcap <- subset(newdata4, Loc2 %in% c('Island'))%>% droplevels()
islcap1<- subset(islcap,islcap$YSW>=	
                   0.00000000 & islcap$YSW <9.191919)
islcap2<- subset(islcap,islcap$YSW>=9.191919)
maincap <- subset(newdata4, Loc2 %in% c('Mainland'))%>% droplevels()


capplot1 <- ggplot() + geom_point(data=rawcap1, aes(YSW,prev, color=Loc2, size=MYSE),alpha=.5, position = "jitter") + geom_ribbon(data = islcap,aes(x=YSW, y=PredictedProb,ymin = LL,ymax = UL), alpha = 0.2) + geom_line(data=islcap1,aes(x=YSW,y=PredictedProb,color = Loc2),linewidth = 1) + geom_line(data=islcap2,aes(x=YSW, y=PredictedProb,colour = Loc2),linewidth = 1, linetype = "dashed") + geom_ribbon(data = maincap,aes(x=YSW, y=PredictedProb,ymin = LL,ymax = UL), alpha = 0.2) + geom_line(data=maincap,aes(x=YSW,y=PredictedProb,color = Loc2),linewidth = 1) + theme_classic() + scale_color_manual(values=c("deepskyblue4","deeppink4"))+ scale_fill_manual(values=c("deepskyblue4","deeppink4")) + theme(text = element_text(size = 20), legend.text=element_text(size=15),legend.position = "top") + scale_size(guide = "none") + theme(legend.title=element_blank()) + xlab("Years since WNS arrival") + ylab("Probability of MYSE capture") + scale_x_continuous(breaks = c(	
  0.0000000, 2.0202020,	4.040404,	6.060606, 8.080808, 10.000000),labels = c("0","2","4","6", "8","10"))
capplot1


