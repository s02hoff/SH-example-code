 
  #
  # Known fate model analysis
  #

  # source the initialization script
  source( here::here( "r", "initscript.r" ) )

  # load data
  duck <- convert.inp( here::here( "data", "BLCKDUCK.inp" ),  ## file names
                         covariates = c("agecov", 
                                        "weight",
                                        "wing",
                                        "condition")         ## names of individual covariates
  )
  
    
  #
  # Process data
  #
  duck.processed <- process.data( duck, model="Known" )
  
  #
  # Create default design data
  #
  duck.ddl <- make.design.data( duck.processed )
  
  
  #
  #  Define models for S
  #
  Sdot <- list( formula=~1 )
  Stime <- list( formula=~time )
  Sage <- list( formula =~ agecov )
  Sweight <- list( formula =~ weight )
  Swing <- list( formula =~ wing )
  Scond <- list( formula =~ condition )
  Swght.con <- list( formula =~ weight + condition )
  Stime.wght <- list( formula =~ weight + time )
  Stime.con <- list( formula =~ time + condition )
  Sage.con <- list( formula =~ agecov + condition )
  
  
  # replace other parameters with assumptions (p=1)
  
  
  
  #
  # Run models
  #
  duck.Sdot <- mark( duck.processed,
                     duck.ddl,
                     model.parameters=list( S=Sdot ),
                     delete = TRUE )
  
  duck.Stime <- mark(duck.processed,
                     duck.ddl,
                     model.parameters=list( S=Stime ),
                     delete = TRUE )
#looks like survival is decreasing over time 
  
  # Now add covariates...  
  
  
  duck.Sage <- mark( duck.processed,
                        duck.ddl,
                        model.parameters=list( S=Sage ),
                        delete = TRUE )
  
  
  duck.Sweight <- mark( duck.processed,
                     duck.ddl,
                     model.parameters=list( S=Sweight ),
                     delete = TRUE )
  
  duck.Swing <- mark( duck.processed,
                        duck.ddl,
                        model.parameters=list( S=Swing ),
                        delete = TRUE )
  
  
  duck.Scond <- mark( duck.processed,
                      duck.ddl,
                      model.parameters=list( S=Scond ),
                      delete = TRUE )
  
  duck.wghtcond <- mark( duck.processed,
                      duck.ddl,
                      model.parameters=list( S=Swght.con ),
                      delete = TRUE )
  
  duck.wghttime <- mark( duck.processed,
                         duck.ddl,
                         model.parameters=list( S=Stime.wght ),
                         delete = TRUE )
  
  duck.contime <- mark( duck.processed,
                         duck.ddl,
                         model.parameters=list( S=Stime.con ),
                         delete = TRUE )
  
  duck.agecon <- mark( duck.processed,
                        duck.ddl,
                        model.parameters=list( S=Sage.con ),
                        delete = TRUE )
  
  model.table()
  
  