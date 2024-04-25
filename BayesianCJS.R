  #
  #Bayesian CJS
  #

  source( here::here( "r", "initscript.r" ) )
  
  # read data file 
  dipper <- read.csv( here::here( "data", "pitNY.csv" ) )
  
  # make encounter history matrix
  EH <- as.matrix( dipper[,1:7] )
  
  nYears <- ncol( EH )
  nAnimal <- nrow( EH )
  
  # determine the occasion of the first capture of each individual
  get.first <- function(x){
  # a function to identify the time period of the first capture from an encounter history matrix
    return( min(which( x != 0 ) ) )
  }
  fall <- apply( EH, 1, get.first )
  
  # remove animals caught only in last year
  EH <- EH[-which(fall==nYears),]
  f <- fall[-which(fall==nYears)]
  
  nAnimal <- nrow( EH )
  
  # Define a list of data to be passed to JAGS
  cjs.data <- list( y=EH, f=f, nind=nAnimal, nocc=nYears )
  
  # Define a function to generate the initial values and return them as a list
  ch.init <- function( ch, f ){
    for( i in 1:dim(ch)[1] ){
      ch[i,1:f[i]] <- NA
    }
    return(ch)
  }
  z <- ch.init( EH, f )    # z is partially observed, so has to agree with y
  z <- ifelse( !is.na(z), 1, z )
  cjs.inits <- function(){
    list(
      z = z,
      b0.phi = runif(1, -3, 3),
      b0.p = runif(1, -3, 3)
    )
  }
  
  # set parameters to track in JAGS
  cjs.parms <- c("b0.phi", "b0.p","mean.phi", "mean.p")
  
  # set up for MCMC run
  ni <- 10000
  nt <- 1
  nb <- 5000
  nc <- 1
  
  # run the MCMC chain in JAGS
  cjs.result <- jags( cjs.data, 
                      cjs.inits,
                      cjs.parms,
                      "jags/cjs_phidot_pdot.txt",
                      n.chains=nc, 
                      n.iter=ni, 
                      n.burnin=nb,
                      n.thin=nt
                    )


