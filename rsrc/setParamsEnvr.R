############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################
# WARNING: requires jsonlite library:
# library(jsonlite)

setParamsEnvr <- function( N.CORES = NULL )
{
  # initialise
  n.cores          <- if( is.null(N.CORES) ) detectCores() else N.CORES
  dt.hst.params <- NA
  lma.params    <- NA
  ini.guess.df  <- NA
  n.guess       <- 0L
  
  # ----------------- internal functions ---------------------------
  .makeDFTableOfInis <- function( ..., Ignore.Seed = FALSE )
  {
    args <- list(...)
    if( !is.null(args[["INI.GUESS.LST"]]) ) ini.guess.params <- args[["INI.GUESS.LST"]]
    else if( !is.null(args[["PATH.TO.INI.GUESS.JSON"]]) )
      ini.guess.params <- read_json(args[["PATH.TO.INI.GUESS.JSON"]], simplifyVector = T)

    n.guess      <<- as.integer( ini.guess.params[["n_guesses"]] )
    ini.guess.df <<- as.data.frame( array( NA, dim = c(n.guess, 12) ) )
    
    names( ini.guess.df ) <<- c("Idir", "IpPs","IoPs1", "IoPs2", 
                                "Sd", "Mu", 
                                "Tau_dir", "Tau_oPs1", "Tau_oPs2", 
                                "Bounded", "PosConstrErr", "UpThrldDt" ) 

    if( !Ignore.Seed ) set.seed( ini.guess.params[["seed"]] )
    # intensities
    int.rng <- ini.guess.params[["i_dir_range_au"]]
    ini.guess.df[["Idir"]] <<- runif( n.guess, int.rng[1], int.rng[2] )    
    ifrac.pps <- 1L
    ifrac.ops <- 1L
    if( ini.guess.params[["preserve_frac_for_ints"]] ){
      ifrac.pps <- 1/6
      ifrac.ops <- 1.5/6 # WARNING: not 3/6 as initially defined for 2 oPs components
    }
    ini.guess.df[["IpPs"]]  <<- ifrac.pps * runif( n.guess, int.rng[1], int.rng[2] )    
    ini.guess.df[["IoPs1"]] <<- ifrac.ops * runif( n.guess, int.rng[1], int.rng[2] )    
    ini.guess.df[["IoPs2"]] <<- ifrac.ops * runif( n.guess, int.rng[1], int.rng[2] )    
    
    # Gaussian part for EMG function
    gauss.rng            <- ini.guess.params[["gauss_range_ns"]]
    ini.guess.df[["Sd"]] <<- runif( n.guess, gauss.rng[1], gauss.rng[2])    
    if( ini.guess.params[["symmetric_mu"]] )
      ini.guess.df[["Mu"]] <<- runif( n.guess, - gauss.rng[2], gauss.rng[2] )    
    else ini.guess.df[["Mu"]] <<- runif( n.guess, gauss.rng[1], gauss.rng[2] )    
    
    # tau's
    ini.guess.df[["Tau_dir"]] <<- runif( n.guess, 
                                         ini.guess.params[["tau_dir_range_ns"]][1], 
                                         ini.guess.params[["tau_dir_range_ns"]][2] )
    ini.guess.df[["Tau_oPs1"]] <<- runif( n.guess, 
                                         ini.guess.params[["tau_ops1_range_ns"]][1], 
                                         ini.guess.params[["tau_ops1_range_ns"]][2] )
    ini.guess.df[["Tau_oPs2"]] <<- runif( n.guess, 
                                          ini.guess.params[["tau_ops2_range_ns"]][1], 
                                          ini.guess.params[["tau_ops2_range_ns"]][2] ) 
    
    # boundaries
    if( ini.guess.params[["no_lower_bounds"]] ) ini.guess.df[["Bounded"]] <<- rep( TRUE, n.guess )
    else ini.guess.df[["Bounded"]] <<- as.logical(sample(c(TRUE, FALSE), n.guess, replace = TRUE))
    
    ini.guess.df[["PosConstrErr"]] <<- 10^rnorm( n.guess, 
                                                 mean = ini.guess.params[["positivity_constr_error"]][1], 
                                                 sd   = ini.guess.params[["positivity_constr_error"]][2] )
    ini.guess.df[["UpThrldDt"]] <<- rlnorm( n.guess,
                                            ini.guess.params[["dt_max_range_for_fast_comp_lrnorm_ns"]][1],
                                            ini.guess.params[["dt_max_range_for_fast_comp_lrnorm_ns"]][2] )
  }  
  
  # -------------- main function ------------------------------------
  # from three files
  setParams <- function( PATH.TO.HST.PARAMS, 
                         PATH.TO.INI.GUESS.PARAMS, 
                         PATH.TO.LMA.PARAMS, 
                         Ignore.Seed = FALSE )
  {
    .makeDFTableOfInis( PATH.TO.INI.GUESS.JSON = PATH.TO.INI.GUESS.PARAMS, Ignore.Seed )
    dt.hst.params <<- read_json( PATH.TO.HST.PARAMS, simplifyVector = T )
    lma.params    <<- read_json( PATH.TO.LMA.PARAMS, simplifyVector = T )
  }
  # from one file
  setParamsCombinedJSON <- function( PATH.TO.PARAMS, 
                                     Ignore.Seed = FALSE )
  {
    params <- read_json( PATH.TO.PARAMS, simplifyVector = T )
    .makeDFTableOfInis( INI.GUESS.LST = params[["ini_guess_params"]], Ignore.Seed )
    dt.hst.params <<- params[["dt_histogram_params"]]
    lma.params    <<- params[["lma_params"]]
  }
  
  rm( list = c("N.CORES") )
  return( environment() )
}