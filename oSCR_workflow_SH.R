  #
  # Spatial capture-recapture analysis
  #

  # if not already installed
  library( devtools )
  install_github( "jaroyle/oSCR" )
  
  # initialization script
  source( here::here( "r", "initscript.r" ) )
  
  # load data
  data( ocelot )

  # examine the data
  ls()
  head(edf1)
  head(edf2)
  head(tdf1)
  head(tdf2)
  
  
  

# Pre-Disturbance Analysis ------------------------------------------------

  
  
  # set up oSCR data frame
  data <- data2oscr( edf = edf1,          # the edf
                     tdf = list(tdf1),    # the tdf
                     sess.col = 1,        # session col NUMBER (edf)
                     id.col = 2,          # ind ID col NUMBER (edf)
                     occ.col = 3,         # occasion col NUMBER (edf)
                     trap.col = 4,        # detector col NUMBER (edf)
                     sex.col = 5,         # sex col NUMBER
                     sex.nacode = "U",    # character for unknown sex?
                     K = 44,              # number of occasions
                     ntraps = 23,         # number of traps
                     remove.extracaps = T)# Encure binary data
  
  # extract the scrFrame
  ocelot.sf <- data$scrFrame
  
  # set buffer and resolution
  ocelot.ss <- make.ssDF( ocelot.sf, 
                          buffer = 2500, 
                          res = 500 )
  
  #can also do make.ssDF in GIS and import raster
  
  #bigger buffer the more the activity centers to predict, want to make it big enough so animals are not up against buffer size
  #any animal within grid has to live within buffer
  #res is tradeoff between res you think is useful given the spacing of the cameras / computer time to run models - smaller resolution means more activity centers, but more computing power
  
  # examine data summary
  ocelot.sf
  
  # plot observations
  plot(ocelot.sf)
  #average coordinates of capture, weighted by number of capture (spatial average of where cat was caught)
  
  # fit model
  scr0.pre <- oSCR.fit( model = list( D~1, p0~1, sig~1 ), 
                        scrFrame = ocelot.sf, 
                        ssDF = ocelot.ss,
                        trimS = 2500)
  
  scr0.pre.psex.sigsex <- oSCR.fit( model = list( D~1, p0~sex, sig~sex  ),
                        scrFrame = ocelot.sf,
                        ssDF = ocelot.ss,
                        trimS = 2500)
  
  scr0.pre.p1.sigsex <- oSCR.fit( model = list( D~1, p0~1, sig~sex  ),
                               scrFrame = ocelot.sf,
                               ssDF = ocelot.ss,
                               trimS = 2500)

  scr0.pre.psex.sig1 <- oSCR.fit( model = list( D~1, p0~sex, sig~1  ),
                                  scrFrame = ocelot.sf,
                                  ssDF = ocelot.ss,
                                  trimS = 2500)  
  
  
  #intercept only model right now, d0 = baseline density (per pixel in this case, set by resolution - 500x500), p0 is prob if animal is right on trap (lamda), and sigma (declining activity around site), psi has to do with sex ratio (sex built into likelihood of these models)
  #trim will help it run faster, can base it off mean max distance 
  
  
  # look at the model summary (shows logit transformed parameters)
  scr0.pre
  scr0.pre.psex.sigsex
  scr0.pre.p1.sigsex
  #p is statistically sig but not strong
  scr0.pre.psex.sig1

  #compile models for viewing model select
  v <- fitList.oSCR(list( scr0.pre, scr0.pre.psex.sigsex, scr0.pre.p1.sigsex,scr0.pre.psex.sig1), rename = T )
  modSel.oSCR(v)
  
  #based on AIC top model is sex difference on sigma (D~1,p0~1,sig~sex)
  #not much better than others due to low sample size for females (not captured on many traps)
  

  
  # estimated density at the pixel scale
  (pre.pix.density <- exp(scr0.pre$outStats$mle[3]))
  (pre.pix.density2 <- exp(scr0.pre.psex.sigsex$outStats$mle[3]))
  (pre.pix.density3 <- exp(scr0.pre.p1.sigsex$outStats$mle[3]))
  (pre.pix.density4 <- exp(scr0.pre.psex.sig1$outStats$mle[3]))
  
  # abundance
  (pre.total.N <- pre.pix.density * nrow(ocelot.ss[[1]])) 
  (pre.total.N2 <- pre.pix.density2 * nrow(ocelot.ss[[1]]))
  (pre.total.N3 <- pre.pix.density3 * nrow(ocelot.ss[[1]]))
  (pre.total.N4 <- pre.pix.density4 * nrow(ocelot.ss[[1]]))
  
  # density
  (pre.total.D <- pre.total.N / (0.5 * 0.5 * nrow(ocelot.ss[[1]])))
  (pre.total.D2 <- pre.total.N2 / (0.5 * 0.5 * nrow(ocelot.ss[[1]])))
  (pre.total.D3 <- pre.total.N3 / (0.5 * 0.5 * nrow(ocelot.ss[[1]])))
  (pre.total.D4 <- pre.total.N4 / (0.5 * 0.5 * nrow(ocelot.ss[[1]])))
  
  
  # density with precision esitmates
  #pixel density (d.factor=1)
  #total density (d.factor=nrow(ocelot.ss[[1]])) #dfactor is multiplier, so use number of rows 
  #also animals per study area / abundance
  #100km density (d.factor=400))
  
  dens.pred.df   <- data.frame(Session = 1)
  (pre.pix.D <- get.real( scr0.pre, type="dens", 
                          newdata = dens.pred.df,
                          d.factor = 1 ) )
  
  (pre.tot.D <- get.real( scr0.pre, type="dens", 
                          newdata = dens.pred.df,
                          d.factor = nrow(ocelot.ss[[1]]) ) )
  
  (pre.km.D  <- get.real( scr0.pre, type="dens", 
                          newdata = dens.pred.df,
                          d.factor = 400 ) )
  
  #density with precision for best model
  
  dens.pred.df   <- data.frame(Session = 1)
  (pre.pix.D <- get.real( scr0.pre.p1.sigsex, type="dens", 
                          newdata = dens.pred.df,
                          d.factor = 1 ) )
  
  (pre.tot.D <- get.real( scr0.pre.p1.sigsex, type="dens", 
                          newdata = dens.pred.df,
                          d.factor = nrow(ocelot.ss[[1]]) ) )
  
  (pre.km.D  <- get.real( scr0.pre.p1.sigsex, type="dens", 
                          newdata = dens.pred.df,
                          d.factor = 400 ) )
  
  #for best model, shows density estimate of 25.10 for females and 47.44 for males
  
  
  
  # baseline detection probability
  (pre.prob.det <- plogis( scr0.pre$outStats$mle[1] ) )
  #daily basis for one camera prob is 0.027
  #
  (pre.prob.det2 <- plogis( scr0.pre.p1.sigsex$outStats$mle[1] ) )
  #detection prob according to best model is 0.03
  
  

  
  
  # sigma
  (pre.sigma <- exp(scr0.pre$outStats$mle[2]))
  #630m, moving big distances (units depends on coorinate system)
  
  sig.pred.df <- data.frame( Session = 1, 
                             sex = c( "M", "F" ) )
  
  (pre.sig <- get.real( scr0.pre, 
                        type="sig", 
                        newdata = sig.pred.df ) )
  
  
  
  
  # sigma for best model
  (pre.sigma <- exp(scr0.pre.p1.sigsex$outStats$mle[2]))
  #699m (units depends on coorinate system)
  
  sig.pred.df2 <- data.frame( Session = 1, 
                             sex = c( "male", "female" ) )
  
  pre.sig2 <- get.real( scr0.pre.p1.sigsex, 
                        type="sig", 
                        newdata = sig.pred.df2 ) 
  
  
  # estimate       se      lwr      upr Session    sex
  # 1 450.5674 76.92780 322.4237 629.6405       1   male
  # 2 699.1881 93.88214 537.4010 909.6819       1 female
  
  

  
  


# Disturbance Analysis ----------------------------------------------------

  
  # set up oSCR data frame
  data1 <- data2oscr( edf = edf2,          # the edf
                     tdf = list(tdf2),    # the tdf
                     sess.col = 1,        # session col NUMBER (edf)
                     id.col = 2,          # ind ID col NUMBER (edf)
                     occ.col = 3,         # occasion col NUMBER (edf)
                     trap.col = 4,        # detector col NUMBER (edf)
                     sex.col = 5,         # sex col NUMBER
                     sex.nacode = "U",    # character for unknown sex?
                     K = 96,              # number of occasions
                     ntraps = 23,         # number of traps
                     remove.extracaps = T)# Encure binary data
  
  
  # extract the scrFrame
  ocelot.sf1 <- data1$scrFrame
  
  # set buffer and resolution
  ocelot.ss1 <- make.ssDF( ocelot.sf1, 
                          buffer = 2500, 
                          res = 500 )
  
  
  
  #run model
  scr0.post.p1.sigsex <- oSCR.fit( model = list( D~1, p0~1, sig~sex  ),
                                  scrFrame = ocelot.sf1,
                                  ssDF = ocelot.ss1,
                                  trimS = 2500)
  
  
  #look at model summary
  scr0.post.p1.sigsex
  
  # estimated density at the pixel scale
  (post.pix.density <- exp(scr0.post.p1.sigsex$outStats$mle[3]))

  
  # abundance
  (post.total.N <- post.pix.density * nrow(ocelot.ss1[[1]])) 

  
  # density
  (post.total.D <- post.total.N / (0.5 * 0.5 * nrow(ocelot.ss[[1]])))

  #density is now 3.12
  
  
  post.dens.pred.df   <- data.frame(Session = 2)
  (post.pix.D <- get.real( scr0.post.p1.sigsex, type="dens", 
                          newdata = post.dens.pred.df,
                          d.factor = 1 ) )
  
  (post.tot.D <- get.real( scr0.post.p1.sigsex, type="dens", 
                          newdata = post.dens.pred.df,
                          d.factor = nrow(ocelot.ss[[1]]) ) )
  
  (post.km.D  <- get.real( scr0.post.p1.sigsex, type="dens", 
                          newdata = post.dens.pred.df,
                          d.factor = 400 ) )
  
  #Pre - animals in study area is 21 and 41, post is 22 and 25
  
  
  
  # baseline detection probability
  (post.prob.det <- plogis( scr0.post.p1.sigsex$outStats$mle[1] ) )
  #daily basis for one camera prob is 0.042

  
  
  
  # sigma
  (post.sigma <- exp(scr0.post.p1.sigsex$outStats$mle[2]))
  #655m
  
  post.sig.pred.df <- data.frame( Session = 2, 
                             sex = c( "male", "female" ) )
  
  (post.sig <- get.real( scr0.post.p1.sigsex, 
                        type="sig", 
                        newdata = post.sig.pred.df ) )
 
  
  #results from pre-disturbance
  
  # estimate       se      lwr      upr Session    sex
  # 1 450.5674 76.92780 322.4237 629.6405       1   male
  # 2 699.1881 93.88214 537.4010 909.6819       1 female
  
  
  #results from post-disturbance
  # estimate       se      lwr      upr Session    sex
  # 1 511.2214 52.51441 417.9930 625.2431       2   male
  # 2 655.2013 57.86728 551.0556 779.0299       2 female 
  
  #males moving more, females moving less?
  
  

  
  
  v1 <- fitList.oSCR(list( scr0.pre.p1.sigsex,scr0.post.p1.sigsex), rename = T )
  modSel.oSCR(v1)
  
  
  
  